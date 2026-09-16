/*
 * Block reader that keeps parsing out of the input critical section.
 * See fast_reader.h.
 */

#include "fast_reader.h"
#include <cstring>
#include <unistd.h>

namespace kraken2 {

static const size_t REFILL_CHUNK = 1 << 20;

FastReader::FastReader()
    : scanned_(0), record_count_(0), loaded_records_(0),
      format_(FORMAT_AUTO_DETECT), fault_count_(0) {
  buf_.reserve(16 << 20);
  nl_.reserve(1 << 17);
  records_.reserve(1 << 15);
}

void FastReader::ScanNewlines() {
  const char *base = buf_.data();
  const char *p = base + scanned_;
  const char *end = base + buf_.size();
  while (p < end) {
    const char *nl = (const char *) memchr(p, '\n', end - p);
    if (! nl) break;
    nl_.push_back((uint32_t) (nl - base));
    p = nl + 1;
  }
  scanned_ = buf_.size();
}

void FastReader::TruncateIndex(size_t keep) {
  size_t n = nl_.size();
  while (n > 0 && (size_t) nl_[n - 1] >= keep)
    n--;
  nl_.resize(n);
  scanned_ = keep;
}

// Reads until `bytes` have been appended or the descriptor is exhausted.  A
// single read() on a pipe returns whatever is available, so loop.
bool FastReader::Fill(int fd, StreamCursor &cur, size_t bytes) {
  size_t base = buf_.size();
  buf_.resize(base + bytes);
  size_t got = 0;
  while (got < bytes) {
    ssize_t n = read(fd, buf_.data() + base + got, bytes - got);
    if (n <= 0)
      break;
    got += (size_t) n;
  }
  buf_.resize(base + got);
  if (got == 0) {
    cur.eof = true;
    return false;
  }
  return true;
}

// The format is the first record marker in the stream.  Blank lines ahead of it
// are skipped, as kseq skips them.
static void DetectFormat(StreamCursor &cur, const std::vector<char> &buf) {
  if (cur.format != FORMAT_AUTO_DETECT)
    return;
  for (size_t i = 0; i < buf.size(); i++) {
    char c = buf[i];
    if (c == '\n' || c == '\r' || c == ' ' || c == '\t')
      continue;
    if (c == '@') cur.format = FORMAT_FASTQ;
    else if (c == '>') cur.format = FORMAT_FASTA;
    else errx(EX_DATAERR, "sequence reader - unrecognized file format");
    return;
  }
}

bool PrimeStream(int fd, StreamCursor &cur) {
  char chunk[1 << 16];
  ssize_t n = read(fd, chunk, sizeof(chunk));
  if (n <= 0) {
    cur.eof = true;
    return false;
  }
  cur.carry.assign(chunk, chunk + n);
  DetectFormat(cur, cur.carry);
  return true;
}

void FastReader::Reset(StreamCursor &cur) {
  records_.clear();
  record_count_ = 0;
  loaded_records_ = 0;
  buf_.clear();
  nl_.clear();
  scanned_ = 0;
  buf_.insert(buf_.end(), cur.carry.begin(), cur.carry.end());
  cur.carry.clear();
  format_ = cur.format;
}

// Length of line `j` with a trailing carriage return dropped, as kseq reads it.
static inline size_t LineLen(const std::vector<char> &buf, size_t start, size_t stop) {
  if (stop > start && buf[stop - 1] == '\r')
    stop--;
  return stop - start;
}

// Appends the offset just past each record that ends within the lines scanned
// so far.  A record is complete once the next header is seen, or once its
// quality has reached the length of its sequence.
void FastReader::CollectRecordEnds(std::vector<size_t> &ends,
                                   RecordScan &scan) const {
  for (; scan.next_line < nl_.size(); scan.next_line++) {
    size_t start = scan.next_line == 0 ? 0 : (size_t) nl_[scan.next_line - 1] + 1;
    size_t stop = (size_t) nl_[scan.next_line];
    size_t len = LineLen(buf_, start, stop);
    char c = len ? buf_[start] : '\0';

    // A header line ends the record whose sequence is being read.
    if (scan.state == 1 && (c == '>' || c == '@')) {
      ends.push_back(start);
      scan.state = 0;
    }
    if (scan.state == 0) {
      if (c == '>' || c == '@') {  // anything before a header is skipped
        scan.state = 1;
        scan.seq_len = scan.quals_len = 0;
      }
      continue;
    }
    if (scan.state == 1) {
      if (c == '+')
        scan.state = 2;
      else
        scan.seq_len += len;
      continue;
    }
    scan.quals_len += len;
    if (scan.quals_len >= scan.seq_len) {
      ends.push_back(stop + 1);
      scan.state = 0;
    }
  }
}

bool FastReader::LoadBlock(int fd, StreamCursor &cur, size_t target_bytes,
                           size_t record_multiple) {
  Reset(cur);
  if (record_multiple < 1)
    record_multiple = 1;

  bool got = Fill(fd, cur, target_bytes);
  if (! got && buf_.empty())
    return false;

  // Keep reading past the target until the block holds `record_multiple` whole
  // records or the input ends, so a record longer than the block is taken whole.
  std::vector<size_t> ends;
  RecordScan scan;
  size_t keep = 0, recs = 0;
  for (;;) {
    // At end of input a final record may lack its trailing newline.
    if (! got && ! buf_.empty() && buf_.back() != '\n')
      buf_.push_back('\n');
    ScanNewlines();
    CollectRecordEnds(ends, scan);
    if (got) {
      recs = ends.size();
      recs -= recs % record_multiple;
      keep = recs ? ends[recs - 1] : 0;
    }
    else {
      // The input ends here, so a record still open is the last one and is
      // kept: Parse reports it the way kseq does.
      recs = ends.size() + (scan.state != 0 ? 1 : 0);
      keep = buf_.size();
    }
    if (keep > 0 || ! got)
      break;
    got = Fill(fd, cur, target_bytes);
  }

  if (keep < buf_.size()) {
    cur.carry.assign(buf_.begin() + keep, buf_.end());
    buf_.resize(keep);
    TruncateIndex(keep);
  }
  loaded_records_ = buf_.empty() ? 0 : recs;
  return ! buf_.empty();
}

bool FastReader::LoadRecords(int fd, StreamCursor &cur, size_t records) {
  Reset(cur);
  if (records == 0)
    return false;

  bool live = true;
  std::vector<size_t> ends;
  RecordScan scan;
  size_t keep = 0, recs = 0;
  ScanNewlines();
  CollectRecordEnds(ends, scan);
  while (ends.size() < records && live) {
    live = Fill(fd, cur, REFILL_CHUNK);
    ScanNewlines();
    CollectRecordEnds(ends, scan);
  }
  if (buf_.empty())
    return false;
  if (! live && buf_.back() != '\n') {
    buf_.push_back('\n');
    ScanNewlines();
    CollectRecordEnds(ends, scan);
  }
  if (ends.size() >= records) {
    keep = ends[records - 1];
    recs = records;
  }
  else {
    keep = buf_.size();
    recs = ends.size() + (scan.state != 0 ? 1 : 0);
  }

  if (keep < buf_.size()) {
    cur.carry.assign(buf_.begin() + keep, buf_.end());
    buf_.resize(keep);
    TruncateIndex(keep);
  }
  loaded_records_ = recs;
  return ! buf_.empty();
}

// A record carries the format kseq would report for it: a quality string makes
// it FASTQ, and its absence makes it FASTA whatever the rest of the file is.
void FastReader::Emit(SeqView &v) {
  v.format = v.quals_len > 0 ? FORMAT_FASTQ : FORMAT_FASTA;
  records_.push_back(v);
}

// A malformed record is still emitted, so the records of a file keep their
// positions and the mates of a pair stay in step.  Its bases are intact, and
// the quality values it does carry still describe the bases they cover, so they
// are kept for quality masking, which bounds itself by the shorter span.  The
// record is written out as FASTA, because a quality string of the wrong length
// would make the classified and unclassified files something Kraken 2 will not
// read back, and padding one to fit would fabricate quality nobody measured.
// The caller reports the fault once the run has finished.
void FastReader::Fault(SeqView &v, const char *verb) {
  fault_count_++;
  if (fault_.empty())
    RecordFault(v, verb);
  Emit(v);
  records_.back().format = FORMAT_FASTA;
}

void FastReader::RecordFault(const SeqView &v, const char *verb) {
  char msg[512];
  snprintf(msg, sizeof(msg),
           "sequence reader - record '%.*s' %s %u bases and %u quality values",
           (int) (v.header_len > 200 ? 200 : v.header_len), v.header, verb,
           v.seq_len, v.quals_len);
  fault_ = msg;
}

void FastReader::Parse() {
  records_.clear();
  char *base = buf_.data();

  auto trim = [](const char *s, const char *e) -> uint32_t {
    if (e > s && e[-1] == '\r')
      e--;
    return (uint32_t) (e - s);
  };
  auto line_start = [&](size_t j) -> char * {
    return base + (j == 0 ? 0 : (size_t) nl_[j - 1] + 1);
  };
  auto line_stop = [&](size_t j) -> char * { return base + nl_[j]; };

  // Splits a header line into identifier and comment the way kseq does: the
  // identifier ends at the first whitespace of any kind.
  auto space = [](char c) {
    return c == ' ' || c == '\t' || c == '\r' || c == '\v' || c == '\f';
  };
  auto set_header = [&](SeqView &v, char *hs, char *he) {
    char *p = hs;
    while (p < he && ! space(*p))
      p++;
    v.header = hs;
    v.header_len = (uint32_t) (p - hs);
    while (p < he && space(*p))
      p++;
    v.comment = p;
    v.comment_len = trim(p, he);
  };

  size_t nlines = nl_.size();
  int state = 0;
  SeqView v;
  char *seq_dst = nullptr, *quals_dst = nullptr;

  for (size_t j = 0; j < nlines; j++) {
    char *ls = line_start(j);
    uint32_t len = trim(ls, line_stop(j));
    char c = len ? *ls : '\0';

    // A header line ends the record being read, which then has no quality, the
    // ordinary case for FASTA and what kseq returns for a FASTQ record cut
    // short before its '+' line.
    if (state == 1 && (c == '>' || c == '@')) {
      v.quals = nullptr;
      v.quals_len = 0;
      Emit(v);
      state = 0;
    }
    if (state == 0) {
      if (c != '>' && c != '@')
        continue;  // kseq skips anything before a header
      v = SeqView();
      set_header(v, ls + 1, ls + len);
      v.seq = v.quals = nullptr;
      v.seq_len = v.quals_len = 0;
      seq_dst = quals_dst = nullptr;
      state = 1;
      continue;
    }
    if (state == 1) {
      if (c == '+') {
        state = 2;
        continue;
      }
      // Sequence spans every line up to the next marker, spliced in place by
      // shifting bytes down over the newlines.
      if (seq_dst == nullptr)
        v.seq = seq_dst = ls;
      if (seq_dst != ls)
        memmove(seq_dst, ls, len);
      seq_dst += len;
      v.seq_len = (uint32_t) (seq_dst - v.seq);
      continue;
    }
    if (quals_dst == nullptr)
      v.quals = quals_dst = ls;  // kseq reads a line here even for an empty read
    if (quals_dst != ls)
      memmove(quals_dst, ls, len);
    quals_dst += len;
    v.quals_len = (uint32_t) (quals_dst - v.quals);
    if (v.quals_len >= v.seq_len) {
      if (v.quals_len != v.seq_len)
        Fault(v, "has");  // kept, so later records keep their positions
      else
        Emit(v);
      state = 0;
    }
  }

  // Whatever is still open belongs to the last record in the input.
  if (state == 1) {
    v.quals = nullptr;
    v.quals_len = 0;
    Emit(v);
  }
  else if (state == 2) {
    if (v.quals_len != v.seq_len)
      Fault(v, "ends with");
    else
      Emit(v);
  }
  record_count_ = records_.size();
}

}  // end namespace
