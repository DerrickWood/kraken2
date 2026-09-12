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
      format_(FORMAT_AUTO_DETECT) {
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

// Appends the offsets of line-initial '>' found since `next_line`, which is a
// line index into the newline index and is advanced past what was examined.  A
// line whose first byte has not been read yet is left for the next call.
void FastReader::CollectHeaders(std::vector<size_t> &heads,
                                size_t &next_line) const {
  for (; next_line <= nl_.size(); next_line++) {
    size_t start = next_line == 0 ? 0 : (size_t) nl_[next_line - 1] + 1;
    if (start >= buf_.size())
      break;
    if (buf_[start] == '>')
      heads.push_back(start);
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
  DetectFormat(cur, buf_);
  format_ = cur.format;

  // Keep reading past the target until the block holds `record_multiple` whole
  // records or the input ends, so a record longer than the block is taken whole.
  std::vector<size_t> heads;
  size_t next_line = 0;
  size_t keep = 0, recs = 0;
  for (;;) {
    // At end of input a final record may lack its trailing newline.
    if (! got && ! buf_.empty() && buf_.back() != '\n')
      buf_.push_back('\n');
    ScanNewlines();
    if (format_ == FORMAT_FASTQ) {
      recs = nl_.size() / 4;
      recs -= recs % record_multiple;
      keep = recs ? (size_t) nl_[recs * 4 - 1] + 1 : 0;
    }
    else {
      CollectHeaders(heads, next_line);
      if (got) {
        // The record under the last header may continue past the buffer.
        recs = heads.empty() ? 0 : heads.size() - 1;
        recs -= recs % record_multiple;
        keep = recs ? heads[recs] : 0;
      }
      else {
        recs = heads.size();
        keep = buf_.size();
      }
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
  // The format decides how records are counted, so settle it before counting.
  if (cur.format == FORMAT_AUTO_DETECT && buf_.empty())
    live = Fill(fd, cur, REFILL_CHUNK);
  DetectFormat(cur, buf_);
  format_ = cur.format;
  ScanNewlines();
  size_t keep = 0, recs = 0;
  if (format_ == FORMAT_FASTQ) {
    const size_t want = records * 4;
    while (nl_.size() < want && live) {
      live = Fill(fd, cur, REFILL_CHUNK);
      ScanNewlines();
    }
    if (buf_.empty())
      return false;
    if (! live && buf_.back() != '\n') {
      buf_.push_back('\n');
      ScanNewlines();
    }
    size_t lines = nl_.size() < want ? (nl_.size() / 4) * 4 : want;
    keep = lines ? (size_t) nl_[lines - 1] + 1 : 0;
    recs = lines / 4;
  }
  else {
    // A FASTA record ends only where the next header begins, and may span any
    // number of lines, so read until the header after the last wanted record is
    // in the buffer, then cut in front of it.
    std::vector<size_t> heads;
    size_t next_line = 0;
    CollectHeaders(heads, next_line);
    while (heads.size() <= records && live) {
      live = Fill(fd, cur, REFILL_CHUNK);
      ScanNewlines();
      CollectHeaders(heads, next_line);
    }
    if (buf_.empty())
      return false;
    if (heads.size() > records) {
      keep = heads[records];
      recs = records;
    }
    else {
      if (buf_.back() != '\n') {
        buf_.push_back('\n');
        ScanNewlines();
      }
      keep = buf_.size();
      recs = heads.size();
    }
  }

  if (keep < buf_.size()) {
    cur.carry.assign(buf_.begin() + keep, buf_.end());
    buf_.resize(keep);
    TruncateIndex(keep);
  }
  loaded_records_ = recs;
  return ! buf_.empty();
}

void FastReader::Parse() {
  records_.clear();
  char *base = buf_.data();

  auto trim = [](const char *s, const char *e) -> uint32_t {
    while (e > s && (e[-1] == '\r' || e[-1] == ' ' || e[-1] == '\t'))
      e--;
    return (uint32_t) (e - s);
  };
  auto line_start = [&](size_t j) -> char * {
    return base + (j == 0 ? 0 : (size_t) nl_[j - 1] + 1);
  };
  auto line_stop = [&](size_t j) -> char * { return base + nl_[j]; };

  // Splits a header line into identifier and comment the way kseq does.
  auto set_header = [&](SeqView &v, char *hs, char *he) {
    char *p = hs;
    while (p < he && *p != ' ' && *p != '\t')
      p++;
    v.header = hs;
    v.header_len = trim(hs, p);
    while (p < he && (*p == ' ' || *p == '\t'))
      p++;
    v.comment = p;
    v.comment_len = trim(p, he);
  };

  if (format_ == FORMAT_FASTQ) {
    size_t recs = nl_.size() / 4;
    records_.resize(recs);
    size_t out = 0;
    for (size_t i = 0; i < recs; i++) {
      size_t j = i * 4;
      char *hs = line_start(j), *he = line_stop(j);
      if (hs >= he || *hs != '@')
        break;
      SeqView &v = records_[out];
      v.format = FORMAT_FASTQ;
      set_header(v, hs + 1, he);
      char *q = line_start(j + 1);
      v.seq = q;
      v.seq_len = trim(q, line_stop(j + 1));
      q = line_start(j + 3);
      v.quals = q;
      v.quals_len = trim(q, line_stop(j + 3));
      out++;
    }
    records_.resize(out);
    record_count_ = out;
    return;
  }

  // FASTA: sequence spans every line up to the next header, spliced in place by
  // shifting bytes down over the newlines.
  size_t nlines = nl_.size();
  size_t j = 0;
  while (j < nlines) {
    char *hs = line_start(j), *he = line_stop(j);
    if (hs >= he || *hs != '>')
      break;
    SeqView v;
    v.format = FORMAT_FASTA;
    v.quals = nullptr;
    v.quals_len = 0;
    set_header(v, hs + 1, he);
    j++;
    char *dst = (j < nlines) ? line_start(j) : base + buf_.size();
    v.seq = dst;
    while (j < nlines) {
      char *ls = line_start(j);
      if (*ls == '>')
        break;
      uint32_t len = trim(ls, line_stop(j));
      if (dst != ls)
        memmove(dst, ls, len);
      dst += len;
      j++;
    }
    v.seq_len = (uint32_t) (dst - v.seq);
    records_.push_back(v);
  }
  record_count_ = records_.size();
}

}  // end namespace
