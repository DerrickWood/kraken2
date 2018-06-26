/*
 * Copyright 2013-2018, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_KRAKEN2_DATA_H_
#define KRAKEN2_KRAKEN2_DATA_H_

#include "kraken2_headers.h"

struct IndexOptions {
  size_t k;
  size_t l;
  uint64_t spaced_seed_mask;
  uint64_t toggle_mask;
  bool dna_db;
};

typedef uint64_t taxid_t;
const taxid_t TAXID_MAX = (taxid_t) -1;

typedef std::unordered_map<taxid_t, uint64_t> taxon_counts_t;

#endif
