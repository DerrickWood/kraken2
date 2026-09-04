/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_AA_TRANSLATE_H_
#define KRAKEN2_AA_TRANSLATE_H_

#include "kraken2_headers.h"

namespace kraken2 {

// Translates a span of bases into all six reading frames.  Takes a pointer and
// length rather than a string so callers holding a view into a read buffer do
// not have to materialize one.
void TranslateToAllFrames(const char *dna_seq, size_t len,
                          std::vector<std::string> &aa_seqs);

}

#endif
