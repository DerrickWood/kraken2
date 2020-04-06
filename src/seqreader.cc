/*
 * Copyright 2013-2020, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include "seqreader.h"

using std::string;

namespace kraken2 {

void StripString(string &str) {
  while (isspace(str.back()))
    str.pop_back();
}

string &Sequence::to_string() {
  str_representation.assign(header);
  switch (format) {
    case FORMAT_FASTQ:
      str_representation.append("\n");
      str_representation.append(seq);
      str_representation.append("\n+\n");
      str_representation.append(quals);
      str_representation.append("\n");
      break;
    default:
      str_representation.append("\n");
      str_representation.append(seq);
      str_representation.append("\n");
      break;
  }
  return str_representation;
}

BatchSequenceReader::BatchSequenceReader() {
  file_format_ = FORMAT_AUTO_DETECT;
  str_buffer_.reserve(8192);
  block_buffer_ = new char[8192];
  block_buffer_size_ = 8192;
}

BatchSequenceReader::~BatchSequenceReader() {
  delete[] block_buffer_;
}

bool BatchSequenceReader::LoadBlock(std::istream &ifs, size_t block_size) {
  ss_.clear();
  ss_.str("");
  if (block_buffer_size_ < block_size) {
    delete[] block_buffer_;
    block_buffer_ = new char[block_size];
    block_buffer_size_ = block_size;
  }
  ifs.read(block_buffer_, block_size);
  if (! ifs && ifs.gcount() <= 0)
    return false;

  if (file_format_ == FORMAT_AUTO_DETECT) {
    switch (block_buffer_[0]) {
      case '@' : file_format_ = FORMAT_FASTQ; break;
      case '>' : file_format_ = FORMAT_FASTA; break;
      default:
        errx(EX_DATAERR, "sequence reader - unrecognized file format");
    }
  }
  str_buffer_.assign(block_buffer_, ifs.gcount());
  ss_ << str_buffer_;
  if (getline(ifs, str_buffer_))
    ss_ << str_buffer_ << "\n";
  if (file_format_ == FORMAT_FASTQ) {
    while (getline(ifs, str_buffer_)) {
      ss_ << str_buffer_ << "\n";
      if (str_buffer_[0] == '@')
        break;
    }
    int lines_to_read = 0;
    if (getline(ifs, str_buffer_)) {
      ss_ << str_buffer_ << "\n";
      lines_to_read = str_buffer_[0] == '@' ? 3 : 2;
      while (lines_to_read-- > 0 && getline(ifs, str_buffer_))
        ss_ << str_buffer_ << "\n";
    }
  }
  else {
    while (ifs) {
      if (ifs.peek() == '>')
        break;
      if (getline(ifs, str_buffer_))
        ss_ << str_buffer_ << "\n";
    }
  }
  return true;
}

bool BatchSequenceReader::LoadBatch(std::istream &ifs, size_t record_count) {
  ss_.clear();
  ss_.str("");
  auto valid = false;
  if (file_format_ == FORMAT_AUTO_DETECT) {
    if (! ifs)
      return false;
    switch (ifs.peek()) {
      case '@' : file_format_ = FORMAT_FASTQ; break;
      case '>' : file_format_ = FORMAT_FASTA; break;
      case EOF : return false;
      default:
        errx(EX_DATAERR, "sequence reader - unrecognized file format");
    }
    valid = true;
  }

  size_t line_count = 0;
  while (record_count > 0 && ifs) {
    if (getline(ifs, str_buffer_))
      line_count++;
    valid = true;
    if (file_format_ == FORMAT_FASTQ) {
      if (line_count % 4 == 0)
        record_count--;
    }
    else {
      if (ifs.peek() == '>')
        record_count--;
    }
    ss_ << str_buffer_ << "\n";
  }

  return valid;
}

bool BatchSequenceReader::NextSequence(Sequence &seq) {
  return BatchSequenceReader::ReadNextSequence
           (ss_, seq, str_buffer_, file_format_);
}

bool BatchSequenceReader::ReadNextSequence(std::istream &is, Sequence &seq,
  std::string &str_buffer, SequenceFormat file_format)
{
  if (! getline(is, str_buffer))
    return false;
  StripString(str_buffer);
  if (file_format == FORMAT_AUTO_DETECT) {
    switch (str_buffer[0]) {
      case '@' : file_format = FORMAT_FASTQ; break;
      case '>' : file_format = FORMAT_FASTA; break;
      default:
        errx(EX_DATAERR, "sequence reader - unrecognized file format");
    }
  }
  seq.format = file_format;
  if (seq.format == FORMAT_FASTQ) {
    if (str_buffer.empty()) // Allow empty line to end file
      return false;
    if (str_buffer[0] != '@')
      errx(EX_DATAERR, "malformed FASTQ file (exp. '@', saw \"%s\"), aborting",
           str_buffer.c_str());
  }
  else if (seq.format == FORMAT_FASTA) {
    if (str_buffer[0] != '>')
      errx(EX_DATAERR, "malformed FASTA file (exp. '>', saw \"%s\"), aborting",
           str_buffer.c_str());
  }
  else
    errx(EX_SOFTWARE, "illegal sequence format encountered in parsing");
  seq.header.assign(str_buffer);
  auto first_whitespace_ch = str_buffer.find_first_of(" \t\r", 1);
  auto substr_len = first_whitespace_ch;
  if (substr_len != std::string::npos)
    substr_len--;
  if (str_buffer.size() > 1)
    seq.id.assign(str_buffer, 1, substr_len);
  else
    return false;

  if (seq.format == FORMAT_FASTQ) {
    if (! getline(is, str_buffer))
      return false;
    StripString(str_buffer);
    seq.seq.assign(str_buffer);
    if (! getline(is, str_buffer))  //  + line, discard
      return false;
    if (! getline(is, str_buffer))
      return false;
    StripString(str_buffer);
    seq.quals.assign(str_buffer);
  }
  else if (seq.format == FORMAT_FASTA) {
    seq.quals.assign("");
    seq.seq.assign("");
    while (is && is.peek() != '>') {
      if (! getline(is, str_buffer))
        return ! seq.seq.empty();
      StripString(str_buffer);
      seq.seq.append(str_buffer);
    }
  }
  return true;
}

}  // end namespace
