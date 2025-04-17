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

  Sequence() : format(SequenceFormat::FORMAT_AUTO_DETECT) {}

  Sequence(Sequence &&other)
      : format(other.format), header(std::move(other.header)),
        comment(std::move(other.comment)), seq(std::move(other.seq)),
        quals(std::move(other.quals))
  {
    // format = other.format;
    // header = std::move(other.header);
    // comment = std::move(other.comment);
    // seq = std::move(other.seq);
    // quals = std::move(other.quals);
  }

  Sequence &operator=(Sequence &other) {
    format = other.format;
    header = std::move(other.header);
    comment = std::move(other.comment);
    seq = std::move(other.seq);
    quals = std::move(other.quals);
    return *this;
  }
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
    size_ = 0;
    // max_size_ = 0;
    file_format_ = SequenceFormat::FORMAT_AUTO_DETECT;
    primary_ = true;
  }

  ~BatchSequenceReader() {
    if (primary_) {
      kseq_destroy(kseq_);
      if (fd_ > 2) {
        close(fd_);
      }
    }
  }

  BatchSequenceReader(const BatchSequenceReader &rhs) {
    kseq_ = rhs.kseq_;
    curr_ = rhs.curr_;
    size_ = rhs.size_;
    // max_size_ = rhs.max_size_;
    // seqs_ = rhs.seqs_;
    primary_ = false;
  }

  bool LoadBlock(size_t block_size) {
    size_t total = 0;
    size_ = 0;
    while (total < block_size) {
      if (kseq_read(kseq_) >= 0) {
        if (size_ == seqs_.size()) {
          seqs_.resize(seqs_.size() + 1);
        }
        copy_from_kseq(seqs_[size_]);
        total += kseq_->seq.l;
        size_++;
      } else {
        break;
      }
    }

    curr_ = 0;
    return size_ > 0;
  }

  bool NextSequence(Sequence& seq)
  {
    bool seq_read = false;
    if (size_ > 0 && curr_ < size_) {
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
    if (size_ > 0 && curr_ < size_) {
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

    // seqs_.resize(i);
    curr_ = 0;
    size_ = i;
    return size_ > 0;
  }

  SequenceFormat file_format() { return file_format_; }

private:
  void copy_from_kseq(Sequence& seq)
  {
    seq.header.assign(kseq_->name.s, kseq_->name.l);
    seq.comment.assign(kseq_->comment.s, kseq_->comment.l);
    seq.seq.assign(kseq_->seq.s, kseq_->seq.l);
    seq.seq.shrink_to_fit();

    if (kseq_->qual.l > 0) {
      file_format_ = SequenceFormat::FORMAT_FASTQ;
      seq.format = SequenceFormat::FORMAT_FASTQ;
      seq.quals.assign(kseq_->qual.s, kseq_->qual.l);
    } else {
      file_format_ = SequenceFormat::FORMAT_FASTA;
      seq.format = SequenceFormat::FORMAT_FASTA;
    }

  }

  kseq_t* kseq_;
  std::vector<Sequence> seqs_;
  SequenceFormat file_format_;
  int fd_;
  size_t curr_;
  size_t size_;
  // size_t max_size_;
  bool primary_;
};

} // namespace kraken2

#endif
