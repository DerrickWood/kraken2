/*
 * Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_KRAKEN2_DATA_H_
#define KRAKEN2_KRAKEN2_DATA_H_

#include "kraken2_headers.h"
#include "readcounts.h"

namespace kraken2 {

struct IndexOptions {
  size_t k;
  size_t l;
  uint64_t spaced_seed_mask;
  uint64_t toggle_mask;
  bool dna_db;
  uint64_t minimum_acceptable_hash_value;
  int revcom_version;   // Fix bug from before K2.0.8
  int db_version;  // To allow for future database structural changes
  int db_type;  // To allow for future use of other data structures
};

typedef uint64_t taxid_t;
const taxid_t TAXID_MAX = (taxid_t) -1;

typedef std::unordered_map<taxid_t, uint64_t> taxon_counts_t;

#ifdef EXACT_COUNTING
  typedef ReadCounts<unordered_set<uint64_t>> READCOUNTER;
#else
  typedef ReadCounts<HyperLogLogPlusMinus<uint64_t>> READCOUNTER;
#endif

typedef std::unordered_map<taxid_t, READCOUNTER> taxon_counters_t;

}

#endif
