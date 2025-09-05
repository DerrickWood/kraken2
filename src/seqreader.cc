/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include "seqreader.h"

using std::string;

namespace kraken2 {

void StripString(string &str) {
  if (str.size() == 0)
    return;
  while (isspace(str.back()))
    str.pop_back();
}

string &Sequence::to_string() {
  str_representation.clear();
  switch (format) {
  case FORMAT_FASTQ:
    str_representation.append("@");
    str_representation.append(header);
    str_representation.append(" ");
    str_representation.append(comment);
    str_representation.append("\n");
    str_representation.append(seq);
    str_representation.append("\n+\n");
    str_representation.append(quals);
    str_representation.append("\n");
    break;
  default:
    str_representation.append(">");
    str_representation.append(header);
    str_representation.append(" ");
    str_representation.append(comment);
    str_representation.append("\n");
    str_representation.append(seq);
    str_representation.append("\n");
    break;
  }
  return str_representation;
}

}  // end namespace
