/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_SEQREADER_H_
#define KRAKEN2_SEQREADER_H_

#include "kraken2_headers.h"
#include "kseq.h"
#include <fcntl.h>

KSEQ_INIT(int, read)

namespace kraken2 {

enum SequenceFormat { FORMAT_AUTO_DETECT,
  FORMAT_FASTA,
  FORMAT_FASTQ };

struct Sequence {
  SequenceFormat format;
  std::string header; // from first char after @/>, up to whitespace
  std::string comment; // from first char after whitespace up to newline
  std::string seq;
  std::string quals; // only meaningful for FASTQ seqs

  std::string& to_string();

  private:
  std::string str_representation;
};

class BatchSequenceReader {
  public:
  BatchSequenceReader(const char* filename = NULL)
  {
    if (filename == NULL) {
      fd_ = fileno(stdin);
    } else {
      fd_ = open(filename, O_RDONLY);
    }

    kseq_ = kseq_init(fd_);
    curr_ = 0;
    file_format_ = SequenceFormat::FORMAT_AUTO_DETECT;
    primary_ = true;
  }

  ~BatchSequenceReader()
  {
    if (primary_) {
      kseq_destroy(kseq_);
      if (fd_ > 2) {
        close(fd_);
      }
    }
  }

  BatchSequenceReader(const BatchSequenceReader& rhs)
  {
    kseq_ = rhs.kseq_;
    curr_ = rhs.curr_;
    seqs_ = rhs.seqs_;
    primary_ = false;
  }

  BatchSequenceReader& operator=(const BatchSequenceReader& rhs) = delete;

  bool LoadBlock(size_t block_size)
  {
    seqs_.resize(0);
    size_t total = 0;
    while (total <= block_size) {
      if (kseq_read(kseq_) >= 0) {
        seqs_.resize(seqs_.size() + 1);
        copy_from_kseq(seqs_.back());

        total += seqs_.back().seq.size();
      } else {
        break;
      }
    }

    curr_ = 0;
    return seqs_.size() > 0;
  }

  bool NextSequence(Sequence& seq)
  {
    bool seq_read = false;
    if (seqs_.size() > 0 && curr_ < seqs_.size()) {
      seq = seqs_[curr_++];
      seq_read = true;
    } else {
      if (kseq_read(kseq_) >= 0) {
        copy_from_kseq(seq);
        seq_read = true;
      }
    }

    return seq_read;
  }

  Sequence* NextSequence()
  {
    if (seqs_.size() > 0 && curr_ < seqs_.size()) {
      return &seqs_[curr_++];
    }

    return NULL;
  }

  bool LoadBatch(size_t record_count)
  {
    size_t i;
    seqs_.resize(record_count);
    for (i = 0; i < record_count; i++) {
      if (kseq_read(kseq_) >= 0) {
        Sequence& seq = seqs_[i];
        copy_from_kseq(seq);
      } else {
        break;
      }
    }

    seqs_.resize(i);
    curr_ = 0;
    return seqs_.size() > 0;
  }

  SequenceFormat file_format() { return file_format_; }

  private:
  void copy_from_kseq(Sequence& seq)
  {
    if (file_format_ == SequenceFormat::FORMAT_AUTO_DETECT) {
      if (kseq_->qual.l > 0)
        file_format_ = SequenceFormat::FORMAT_FASTQ;
      else
        file_format_ = SequenceFormat::FORMAT_FASTA;
    }

    seq.header.assign(kseq_->name.s, kseq_->name.l);
    seq.comment.assign(kseq_->comment.s, kseq_->comment.l);
    seq.seq.assign(kseq_->seq.s, kseq_->seq.l);
  }

  kseq_t* kseq_;
  std::vector<Sequence> seqs_;
  SequenceFormat file_format_;
  int fd_;
  size_t curr_;
  bool primary_;
};

} // namespace kraken2

#endif
