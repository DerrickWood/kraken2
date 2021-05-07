/*
 * Copyright 2013-2020, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_MMSCANNER_H_
#define KRAKEN2_MMSCANNER_H_

#include "kraken2_headers.h"

namespace kraken2 {

const uint64_t DEFAULT_TOGGLE_MASK = 0xe37e28c4271b5a2dULL;
const uint64_t DEFAULT_SPACED_SEED_MASK = 0;
const int BITS_PER_CHAR_DNA = 2;
const int BITS_PER_CHAR_PRO = 4;
const int CURRENT_REVCOM_VERSION = 1;

struct MinimizerData {
  uint64_t candidate;
  ssize_t pos;
};

class MinimizerScanner {
  public:
  // Create scanner for seq over interval [start, finish)
  // Will report length l minimizers for all k-mers in the interval
  // Minimizers ordered by XORing candidates w/ toggle_mask
  //   toggle_mask == 0 implies lexicographical ordering
  MinimizerScanner(ssize_t k, ssize_t l,
                   uint64_t spaced_seed_mask = DEFAULT_SPACED_SEED_MASK,
                   bool dna_sequence = true,
                   uint64_t toggle_mask = DEFAULT_TOGGLE_MASK,
                   int revcom_version = CURRENT_REVCOM_VERSION);

  void LoadSequence(const std::string &seq, size_t start = 0,
      size_t finish = SIZE_MAX);

  uint64_t *NextMinimizer();
  // Return last minimizer, only valid if NextMinimizer last returned non-NULL
  uint64_t last_minimizer() const { return last_minimizer_; }
  ssize_t k() const { return k_; }
  ssize_t l() const { return l_; }
  bool is_dna() const { return dna_; }
  bool is_ambiguous() const {
    return (queue_pos_ < k_ - l_) || (!! last_ambig_);
  }

  private:
  uint64_t reverse_complement(uint64_t kmer, uint8_t n);
  uint64_t canonical_representation(uint64_t kmer, uint8_t n);
  void set_lookup_table_character(char ch, uint8_t val);

  const std::string *str_;  // pointer to sequence
  ssize_t k_;
  ssize_t l_;
  size_t str_pos_, start_, finish_;
  uint64_t spaced_seed_mask_;
  bool dna_;
  uint64_t toggle_mask_;
  uint64_t lmer_;
  uint64_t lmer_mask_;
  uint64_t last_minimizer_;
  ssize_t loaded_ch_;
  std::vector<MinimizerData> queue_;
  ssize_t queue_pos_;
  uint64_t last_ambig_;
  uint8_t lookup_table_[UINT8_MAX + 1];
  const int revcom_version_;
};

}

#endif
