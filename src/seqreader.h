/*
 * Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_SEQREADER_H_
#define KRAKEN2_SEQREADER_H_

#include "kraken2_headers.h"

namespace kraken2 {

enum SequenceFormat {
  FORMAT_AUTO_DETECT,
  FORMAT_FASTA,
  FORMAT_FASTQ
};

struct Sequence {
  SequenceFormat format;
  std::string header;  // header line, including @/>, but not newline
  std::string id;      // from first char. after @/> up to first whitespace
  std::string seq;
  std::string quals;   // only meaningful for FASTQ seqs

  std::string &to_string();

  private:
  std::string str_representation;
};

class BatchSequenceReader {
  public:
  BatchSequenceReader();
  ~BatchSequenceReader();
  BatchSequenceReader(const BatchSequenceReader &rhs) = delete;
  BatchSequenceReader& operator=(const BatchSequenceReader &rhs) = delete;

  bool LoadBatch(std::istream &ifs, size_t record_count);
  bool LoadBlock(std::istream &ifs, size_t block_size);
  bool NextSequence(Sequence &seq);
  static bool ReadNextSequence(std::istream &is, Sequence &seq, 
    std::string &str_buffer_ptr, SequenceFormat format = FORMAT_AUTO_DETECT);

  SequenceFormat file_format() { return file_format_; }

  private:
  std::stringstream ss_;
  std::string str_buffer_;  // used to prevent realloc upon every load/parse
  SequenceFormat file_format_;
  char *block_buffer_;
  size_t block_buffer_size_;
};

}

#endif
