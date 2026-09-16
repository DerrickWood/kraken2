/*
 * Block reader that keeps parsing out of the input critical section.
 * See fast_reader.h.
 */

#include "fast_reader.h"
#include <cerrno>
#include <cstring>
#include <fcntl.h>
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

// Offsets into the buffer and record lengths are 32-bit, so no buffer may grow
// past this.  Only a single record of about 4 GiB can get there.
static const size_t MAX_BUFFER_BYTES = 0xffffffffu;

// One read(), retried when a signal interrupts it, and switched to blocking mode
// once if the descriptor was left non-blocking, so neither is mistaken for the
// end of the input.  Returns what read() returns otherwise.
static ssize_t ReadSome(int fd, char *dst, size_t n) {
  bool made_blocking = false;
  for (;;) {
    ssize_t got = read(fd, dst, n);
    if (got >= 0)
      return got;
    if (errno == EINTR)
      continue;
    if ((errno == EAGAIN || errno == EWOULDBLOCK) && ! made_blocking) {
      int flags = fcntl(fd, F_GETFL);
      if (flags >= 0 && fcntl(fd, F_SETFL, flags & ~O_NONBLOCK) == 0) {
        made_blocking = true;
        continue;
      }
    }
    return got;
  }
}

// Reads until `bytes` have been appended or the descriptor is exhausted.  A
// single read() on a pipe returns whatever is available, so loop.
bool FastReader::Fill(int fd, StreamCursor &cur, size_t bytes) {
  size_t base = buf_.size();
  if (base + bytes > MAX_BUFFER_BYTES) {
    cur.error = "sequence reader - a record is longer than the 4 GiB this "
                "reader supports";
    cur.eof = true;
    return false;
  }
  buf_.resize(base + bytes);
  size_t got = 0;
  while (got < bytes) {
    ssize_t n = ReadSome(fd, buf_.data() + base + got, bytes - got);
    if (n < 0) {
      cur.error = std::string("sequence reader - read error: ") + strerror(errno);
      break;
    }
    if (n == 0)
      break;
    got += (size_t) n;
  }
  buf_.resize(base + got);
  if (got == 0 || ! cur.error.empty()) {
    cur.eof = true;
    return false;
  }
  return true;
}

// Names the compression format whose magic number opens the stream, if any.
// Compressed bytes can contain a newline followed by '@' or '>' by chance, so
// they are recognized by their magic number rather than by the absence of a
// record marker.
static const char *CompressionName(const unsigned char *b, size_t n) {
  if (n >= 2 && b[0] == 0x1f && b[1] == 0x8b) return "gzip";
  if (n >= 3 && b[0] == 'B' && b[1] == 'Z' && b[2] == 'h') return "bzip2";
  if (n >= 6 && b[0] == 0xfd && b[1] == '7' && b[2] == 'z' && b[3] == 'X' &&
      b[4] == 'Z' && b[5] == 0) return "xz";
  if (n >= 4 && b[0] == 0x28 && b[1] == 0xb5 && b[2] == 0x2f && b[3] == 0xfd)
    return "zstd";
  return nullptr;
}

// The format is that of the first line beginning with a record marker.  Only
// blank lines and comment lines, beginning with '#' or ';', may come ahead of
// it; they are skipped, as kseq and the parser here skip them.  Anything else
// first means the input is not sequence data, which is also what keeps a
// document with a line that happens to begin with '>' from being read as FASTA.
static void DetectFormat(StreamCursor &cur, const std::vector<char> &buf) {
  if (cur.format != FORMAT_AUTO_DETECT)
    return;
  for (size_t i = 0; i < buf.size(); ) {
    const char *nl = (const char *) memchr(buf.data() + i, '\n', buf.size() - i);
    size_t stop = nl ? (size_t) (nl - buf.data()) : buf.size();
    size_t j = i;
    while (j < stop && (buf[j] == ' ' || buf[j] == '\t' || buf[j] == '\r'))
      j++;
    if (j < stop) {
      char c = buf[i];
      if (c == '@') { cur.format = FORMAT_FASTQ; return; }
      if (c == '>') { cur.format = FORMAT_FASTA; return; }
      if (c != '#' && c != ';')
        errx(EX_DATAERR, "sequence reader - unrecognized file format");
      for (size_t k = i; k < stop; k++)
        if (buf[k] == '\0')
          errx(EX_DATAERR, "sequence reader - unrecognized file format");
    }
    if (! nl)
      break;
    i = stop + 1;
  }
}

bool PrimeStream(int fd, StreamCursor &cur) {
  char chunk[1 << 16];
  ssize_t n = ReadSome(fd, chunk, sizeof(chunk));
  if (n < 0) {
    cur.error = std::string("sequence reader - read error: ") + strerror(errno);
    cur.eof = true;
    return false;
  }
  if (n == 0) {
    cur.eof = true;
    return false;
  }
  const char *compressed = CompressionName((const unsigned char *) chunk, n);
  if (compressed)
    errx(EX_DATAERR, "sequence reader - input is %s-compressed; decompress it "
                     "before classifying", compressed);
  // A UTF-8 byte order mark, as some editors write, is not part of the data.
  size_t skip = (n >= 3 && (unsigned char) chunk[0] == 0xef &&
                 (unsigned char) chunk[1] == 0xbb &&
                 (unsigned char) chunk[2] == 0xbf) ? 3 : 0;
  cur.carry.assign(chunk + skip, chunk + n);
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
  if (! cur.error.empty()) {
    buf_.clear();
    return false;
  }
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
    if (! got && cur.error.empty() && ! buf_.empty() && buf_.back() != '\n')
      buf_.push_back('\n');
    ScanNewlines();
    CollectRecordEnds(ends, scan);
    if (got || ! cur.error.empty()) {
      // A stream that failed has not ended, so a record still open was cut
      // short by the failure and is left out.
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
    // A failed stream is never read again, so a refused tail, which may be a
    // record of several GiB, is dropped rather than carried.
    if (cur.error.empty())
      cur.carry.assign(buf_.begin() + keep, buf_.end());
    buf_.resize(keep);
    TruncateIndex(keep);
  }
  loaded_records_ = buf_.empty() ? 0 : recs;
  return ! buf_.empty();
}

bool FastReader::LoadRecords(int fd, StreamCursor &cur, size_t records) {
  Reset(cur);
  if (! cur.error.empty()) {
    buf_.clear();
    return false;
  }
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
  if (! live && cur.error.empty() && buf_.back() != '\n') {
    buf_.push_back('\n');
    ScanNewlines();
    CollectRecordEnds(ends, scan);
  }
  if (ends.size() >= records) {
    keep = ends[records - 1];
    recs = records;
  }
  else if (! cur.error.empty()) {
    // Cut short by a failed stream rather than by its end.
    recs = ends.size();
    keep = recs ? ends[recs - 1] : 0;
  }
  else {
    keep = buf_.size();
    recs = ends.size() + (scan.state != 0 ? 1 : 0);
  }

  if (keep < buf_.size()) {
    if (cur.error.empty())
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
    if (p < he)
      p++;  // kseq consumes the one separator and keeps any further whitespace
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
