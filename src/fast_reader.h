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
 * Record boundaries come from counting newlines rather than scanning for a line
 * that begins with '@', which is unsound because Illumina quality strings
 * contain '@' at Phred 31.  Each block is cut on a boundary, so the next begins
 * on one and the count is exact by induction.  Bytes of a trailing partial
 * record are carried forward in a StreamCursor shared by the threads reading
 * that file, so a block costs one read and one small memcpy.
 */

#ifndef KRAKEN2_FAST_READER_H_
#define KRAKEN2_FAST_READER_H_

#include "kraken2_headers.h"
#include "seqreader.h"

namespace kraken2 {

// Per-stream state: the tail of an incomplete record, carried to the next
// block.  One per input file, shared by all threads, only ever touched with the
// input lock held.
struct StreamCursor {
  std::vector<char> carry;
  SequenceFormat format;
  bool eof;
  StreamCursor() : format(FORMAT_AUTO_DETECT), eof(false) { }
};

class FastReader {
  public:
  FastReader();

  // Both loaders must be called with the input lock held, and both leave the
  // reader holding a whole number of records.

  // Pulls about `target_bytes`, then trims back to the last complete record.
  // `record_multiple` forces the kept record count to be a multiple of that
  // value, so interleaved mate pairs are never split across blocks.
  bool LoadBlock(int fd, StreamCursor &cur, size_t target_bytes,
                 size_t record_multiple = 1);

  // Pulls exactly `records` records, for the second mate of a pair so the two
  // files stay in step.
  bool LoadRecords(int fd, StreamCursor &cur, size_t records);

  // Number of records held, without parsing them.
  size_t RecordCount() const;

  // Call outside the lock.  Splits the loaded bytes into records in place.
  void Parse();

  size_t size() const { return records_.size(); }
  const Sequence &operator[](size_t i) const { return records_[i]; }
  Sequence &at(size_t i) { return records_[i]; }
  SequenceFormat file_format() const { return format_; }

  private:
  void ScanNewlines();
  void TruncateIndex(size_t keep);
  bool Fill(int fd, StreamCursor &cur, size_t bytes);

  std::vector<char> buf_;
  std::vector<uint32_t> nl_;   // offsets of newlines within buf_
  size_t scanned_;             // bytes of buf_ already covered by nl_
  std::vector<Sequence> records_;
  size_t record_count_;
  SequenceFormat format_;
};

}  // end namespace

#endif  // KRAKEN2_FAST_READER_H_
