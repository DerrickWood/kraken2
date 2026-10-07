/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_UTILITIES_H_
#define KRAKEN2_UTILITIES_H_

#include "kraken2_headers.h"
#include <vector>
#include <inttypes.h>

// Functions used by 2+ programs that I couldn't think of a better place for.

namespace kraken2 {

// Turns a simple bitstring of 0s and 1s into one where the 1s and 0s are expanded,
// e.g. 010110 expanded by a factor of 2 would become 001100111100
// Allows specification of the spaced seed to represent positions in the sequence
// rather than actual bits in the internal representation
void ExpandSpacedSeedMask(uint64_t &spaced_seed_mask, const int bit_expansion_factor);

std::vector<std::string> SplitString(const std::string &str,
  const std::string &delim = "\t", const size_t max_fields = (size_t) -1);

class BitVec {
public:
  BitVec() : num_ones(0) {}

  void set_bit(size_t pos) {
    size_t num_words = round_to_nearest_word(pos);
    if (num_words > vec.size()) {
      vec.resize(num_words, 0);
    }
    size_t bit = 1;
    vec[pos / block_size] |= bit << (pos % block_size);

    num_ones += 1;
  }

  void merge(const BitVec &other) {
    num_ones = 0;
    if (other.vec.size() > vec.size()) {
      vec.resize(other.vec.size(), 0);
    }
    for (size_t i = 0; i < vec.size(); i++) {
      vec[i] |= other.vec[i];
      num_ones += __builtin_popcountll(vec[i]);
    }
  }

  size_t num_bytes() { return vec.size() * sizeof(size_t); }
  size_t count_ones() { return num_ones; }

  size_t *data() { return vec.data(); }

private:
  size_t round_to_nearest_word(size_t pos) {
    return (pos + block_size + 1) / block_size;
  }

  size_t num_ones;
  std::vector<size_t> vec;
  const size_t block_size = sizeof(size_t) * 8;
};

} // namespace kraken2

#endif
