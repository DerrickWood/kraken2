/*
 * Block reader that keeps parsing out of the input critical section.
 *
 * BatchSequenceReader pulls records with kseq_read and copies each one into a
 * Sequence while the caller holds critical(seqread), so every thread's parsing
 * is serialized against every other thread's.  A sampling profile of classify at
 * 8 threads puts about half of all samples in __psynch_mutexwait.
 *
 * Here the critical section does two things only: read raw bytes, and cut the
 * buffer at an exact record boundary.  Parsing runs afterwards, outside the
 * lock, over the thread's own buffer.  Langmead et al. 2019 (Bioinformatics
 * 35:421) measures this arrangement, which it calls B/L-parsing, as the best of
 * four for thread scaling; kseq inside the lock is their O-parsing.
 *
 * Record boundaries follow the same rules kseq applies, for both formats.  A
 * record begins at a line starting with '>' or '@'.  Its sequence runs over any
 * number of lines and ends at the next line starting with '>', '@' or '+'.  A
 * '+' line begins a quality string, which takes as many lines as it needs to
 * reach the length of the sequence, so a '@' beginning a quality line is never
 * mistaken for a header.  Bytes of a trailing partial record are carried
 * forward in a StreamCursor shared by the threads reading that file, so a block
 * costs one read and one small memcpy.
 */

#ifndef KRAKEN2_FAST_READER_H_
#define KRAKEN2_FAST_READER_H_

#include "kraken2_headers.h"
#include "seqreader.h"

namespace kraken2 {

// A record as pointers into the reader's own buffer, valid until the next load.
// Avoids the four std::string copies a Sequence costs per record; the bases are
// also exactly the flat span the minimizer scanner and a device kernel want.
struct SeqView {
  const char *header;   // identifier, up to the first whitespace
  const char *comment;  // remainder of the header line, may be null
  const char *seq;
  const char *quals;    // null for FASTA
  uint32_t header_len;
  uint32_t comment_len;
  uint32_t seq_len;
  uint32_t quals_len;
  SequenceFormat format;

  // Quality masking rewrites bases in place in the reader's buffer.
  char *seq_mutable() const { return const_cast<char *>(seq); }
};

// Per-stream state: the tail of an incomplete record, carried to the next
// block.  One per input file, shared by all threads, only ever touched with the
// input lock held.
struct StreamCursor {
  std::vector<char> carry;
  SequenceFormat format;
  bool eof;
  StreamCursor() : format(FORMAT_AUTO_DETECT), eof(false) { }
};

// Reads the first bytes of a stream into the cursor's carry and settles the
// stream's format from them, before any thread takes the input lock.  Nothing is
// consumed, since the loaders start from the carry.  Returns false for an empty
// stream, whose format stays FORMAT_AUTO_DETECT.
bool PrimeStream(int fd, StreamCursor &cur);

class FastReader {
  public:
  FastReader();

  // Both loaders must be called with the input lock held, on a cursor that has
  // been through PrimeStream, and both leave the reader holding a whole number
  // of records.

  // Pulls about `target_bytes`, then trims back to the last complete record.
  // `record_multiple` forces the kept record count to be a multiple of that
  // value, so interleaved mate pairs are never split across blocks.  A record
  // longer than `target_bytes` grows the block rather than ending the input.
  bool LoadBlock(int fd, StreamCursor &cur, size_t target_bytes,
                 size_t record_multiple = 1);

  // Pulls exactly `records` records, for the second mate of a pair so the two
  // files stay in step.
  bool LoadRecords(int fd, StreamCursor &cur, size_t records);

  // Number of whole records the last load kept, without parsing them.
  size_t RecordCount() const { return loaded_records_; }

  // Call outside the lock.  Splits the loaded bytes into records in place.
  void Parse();

  // Empty unless a record's quality string and its sequence disagreed in
  // length.  Such a record is still emitted and classified, so the run finishes
  // and the caller reports this at the end.  The count covers every such record
  // this reader has seen, not only the first.
  const std::string &fault() const { return fault_; }
  size_t fault_count() const { return fault_count_; }

  size_t size() const { return records_.size(); }
  const SeqView &operator[](size_t i) const { return records_[i]; }
  SeqView &at(size_t i) { return records_[i]; }
  SequenceFormat file_format() const { return format_; }

  private:
  void ScanNewlines();
  void TruncateIndex(size_t keep);
  bool Fill(int fd, StreamCursor &cur, size_t bytes);
  // Walks lines with kseq's record rules, appending the offset just past each
  // complete record.  State persists across calls so a buffer that grows is
  // scanned once.
  struct RecordScan {
    size_t next_line;  // first line not yet examined
    int state;         // 0 seeking a header, 1 in sequence, 2 in quality
    size_t seq_len, quals_len;
    RecordScan() : next_line(0), state(0), seq_len(0), quals_len(0) { }
  };
  void CollectRecordEnds(std::vector<size_t> &ends, RecordScan &scan) const;

  void Emit(SeqView &v);
  void Fault(SeqView &v, const char *verb);
  void RecordFault(const SeqView &v, const char *verb);

  void Reset(StreamCursor &cur);

  std::vector<char> buf_;
  std::vector<uint32_t> nl_;   // offsets of newlines within buf_
  size_t scanned_;             // bytes of buf_ already covered by nl_
  std::vector<SeqView> records_;
  size_t record_count_;
  size_t loaded_records_;      // whole records kept by the last load
  SequenceFormat format_;
  std::string fault_;
  size_t fault_count_;
};

}  // end namespace

#endif  // KRAKEN2_FAST_READER_H_
