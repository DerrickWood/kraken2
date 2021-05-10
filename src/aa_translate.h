/*
 * Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_AA_TRANSLATE_H_
#define KRAKEN2_AA_TRANSLATE_H_

#include "kraken2_headers.h"

namespace kraken2 {

void TranslateToAllFrames(std::string &dna_seq, std::vector<std::string> &aa_seqs);

}

#endif
