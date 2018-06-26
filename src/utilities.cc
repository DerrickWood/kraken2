/*
 * Copyright 2013-2018, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include "utilities.h"

namespace kraken2 {

void ExpandSpacedSeedMask(uint64_t &spaced_seed_mask, const int bit_expansion_factor) {
  uint64_t new_mask = 0;
  uint64_t bits = 1 << bit_expansion_factor;
  bits--;

  for (int i = 64 / bit_expansion_factor - 1; i >= 0; i--) {
    new_mask <<= bit_expansion_factor;
    if ((spaced_seed_mask >> i) & 1)
      new_mask |= bits;
  }
  spaced_seed_mask = new_mask;
}

}
