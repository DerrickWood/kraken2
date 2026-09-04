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
    : scanned_(0), record_count_(0), format_(FORMAT_AUTO_DETECT) {
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

static void DetectFormat(StreamCursor &cur, const std::vector<char> &buf) {
  if (cur.format != FORMAT_AUTO_DETECT || buf.empty())
    return;
  if (buf[0] == '@') cur.format = FORMAT_FASTQ;
  else if (buf[0] == '>') cur.format = FORMAT_FASTA;
  else errx(EX_DATAERR, "sequence reader - unrecognized file format");
}

bool FastReader::LoadBlock(int fd, StreamCursor &cur, size_t target_bytes,
                           size_t record_multiple) {
  records_.clear();
  record_count_ = 0;
  buf_.clear();
  nl_.clear();
  scanned_ = 0;
  buf_.insert(buf_.end(), cur.carry.begin(), cur.carry.end());
  cur.carry.clear();

  bool got = Fill(fd, cur, target_bytes);
  if (! got && buf_.empty())
    return false;
  DetectFormat(cur, buf_);
  format_ = cur.format;

  // At end of input a final record may lack its trailing newline.
  if (! got && ! buf_.empty() && buf_.back() != '\n')
    buf_.push_back('\n');

  if (record_multiple < 1)
    record_multiple = 1;
  ScanNewlines();
  while (record_multiple > 1 && got && nl_.size() / 4 < record_multiple) {
    got = Fill(fd, cur, target_bytes);
    if (! got && ! buf_.empty() && buf_.back() != '\n')
      buf_.push_back('\n');
    ScanNewlines();
  }

  size_t keep;
  if (format_ == FORMAT_FASTQ) {
    size_t recs = nl_.size() / 4;
    recs -= recs % record_multiple;
    keep = recs ? (size_t) nl_[recs * 4 - 1] + 1 : 0;
  }
  else {
    // FASTA: cut before the last line-initial '>', whose record may continue
    // into the next block.
    keep = buf_.size();
    if (got) {
      size_t i = buf_.size();
      while (i > 0) {
        i--;
        if (buf_[i] == '>' && (i == 0 || buf_[i - 1] == '\n')) {
          keep = i;
          break;
        }
      }
    }
  }

  if (keep < buf_.size()) {
    cur.carry.assign(buf_.begin() + keep, buf_.end());
    buf_.resize(keep);
    TruncateIndex(keep);
  }
  return ! buf_.empty();
}

bool FastReader::LoadRecords(int fd, StreamCursor &cur, size_t records) {
  records_.clear();
  record_count_ = 0;
  buf_.clear();
  nl_.clear();
  scanned_ = 0;
  buf_.insert(buf_.end(), cur.carry.begin(), cur.carry.end());
  cur.carry.clear();
  if (records == 0)
    return false;

  DetectFormat(cur, buf_);
  bool live = true;
  ScanNewlines();
  const size_t want_lines =
      (cur.format == FORMAT_FASTA) ? records * 2 : records * 4;
  while (nl_.size() < want_lines && live) {
    live = Fill(fd, cur, REFILL_CHUNK);
    ScanNewlines();
  }
  if (buf_.empty())
    return false;
  DetectFormat(cur, buf_);
  format_ = cur.format;

  if (! live && ! buf_.empty() && buf_.back() != '\n') {
    buf_.push_back('\n');
    ScanNewlines();
  }

  size_t keep;
  if (format_ == FORMAT_FASTQ) {
    size_t want = records * 4;
    size_t lines = nl_.size() < want ? (nl_.size() / 4) * 4 : want;
    keep = lines ? (size_t) nl_[lines - 1] + 1 : 0;
  }
  else {
    keep = buf_.size();
  }
  if (keep < buf_.size()) {
    cur.carry.assign(buf_.begin() + keep, buf_.end());
    buf_.resize(keep);
    TruncateIndex(keep);
  }
  return ! buf_.empty();
}

size_t FastReader::RecordCount() const {
  if (format_ == FORMAT_FASTQ)
    return nl_.size() / 4;
  size_t n = 0;
  for (size_t i = 0; i < buf_.size(); i++)
    if (buf_[i] == '>' && (i == 0 || buf_[i - 1] == '\n'))
      n++;
  return n;
}

void FastReader::Parse() {
  char *base = buf_.data();
  size_t out = 0;

  auto trim = [](const char *s, const char *e) -> size_t {
    while (e > s && (e[-1] == '\r' || e[-1] == ' ' || e[-1] == '\t'))
      e--;
    return (size_t) (e - s);
  };
  auto line_start = [&](size_t j) -> char * {
    return base + (j == 0 ? 0 : (size_t) nl_[j - 1] + 1);
  };
  auto line_stop = [&](size_t j) -> char * { return base + nl_[j]; };

  // Splits a header line into the identifier and the comment that follows it,
  // matching how kseq fills Sequence::header and Sequence::comment.
  auto set_header = [&](Sequence &s, char *hs, char *he) {
    char *p = hs;
    while (p < he && *p != ' ' && *p != '\t')
      p++;
    s.header.assign(hs, trim(hs, p));
    if (p < he) {
      while (p < he && (*p == ' ' || *p == '\t'))
        p++;
      s.comment.assign(p, trim(p, he));
    }
    else {
      s.comment.clear();
    }
  };

  if (format_ == FORMAT_FASTQ) {
    size_t recs = nl_.size() / 4;
    if (records_.size() < recs)
      records_.resize(recs);
    for (size_t i = 0; i < recs; i++) {
      size_t j = i * 4;
      char *hs = line_start(j), *he = line_stop(j);
      if (hs >= he || *hs != '@')
        break;
      Sequence &s = records_[out];
      s.format = FORMAT_FASTQ;
      set_header(s, hs + 1, he);
      char *q = line_start(j + 1);
      s.seq.assign(q, trim(q, line_stop(j + 1)));
      q = line_start(j + 3);
      s.quals.assign(q, trim(q, line_stop(j + 3)));
      out++;
    }
    record_count_ = out;
    records_.resize(out);
    return;
  }

  // FASTA: a record's sequence spans every line up to the next header.
  size_t nlines = nl_.size();
  size_t j = 0;
  while (j < nlines) {
    char *hs = line_start(j), *he = line_stop(j);
    if (hs >= he || *hs != '>')
      break;
    if (records_.size() <= out)
      records_.resize(out + 1);
    Sequence &s = records_[out];
    s.format = FORMAT_FASTA;
    set_header(s, hs + 1, he);
    s.seq.clear();
    s.quals.clear();
    j++;
    while (j < nlines) {
      char *ls = line_start(j);
      if (*ls == '>')
        break;
      s.seq.append(ls, trim(ls, line_stop(j)));
      j++;
    }
    out++;
  }
  record_count_ = out;
  records_.resize(out);
}

}  // end namespace
