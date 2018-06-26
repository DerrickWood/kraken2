/*
 * Copyright 2013-2018, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include "mmscanner.h"
#include <cassert>

using std::string;
using std::vector;

namespace kraken2 {

#define BITS_PER_CHAR_DNA 2
#define BITS_PER_CHAR_PRO 4

MinimizerScanner::MinimizerScanner(ssize_t k, ssize_t l,
    uint64_t toggle_mask, bool dna_sequence, uint64_t spaced_seed_mask)
    : str_(nullptr), k_(k), l_(l), str_pos_(0), start_(0), finish_(0),
      toggle_mask_(toggle_mask), dna_(dna_sequence),
      spaced_seed_mask_(spaced_seed_mask), loaded_ch_(0),
      last_ambig_(0)
{
  if (l_ > (ssize_t) ((sizeof(uint64_t) * 8 - 1) / (dna_ ? BITS_PER_CHAR_DNA : BITS_PER_CHAR_PRO)))
    errx(EX_SOFTWARE, "l exceeds size limits for minimizer %s scanner",
         dna_ ? "nucleotide" : "protein");
  lmer_mask_ = 1;
  lmer_mask_ <<= (l_ * (dna_ ? BITS_PER_CHAR_DNA : BITS_PER_CHAR_PRO));
  lmer_mask_--;
  toggle_mask_ &= lmer_mask_;
  if (finish_ == SIZE_MAX)
    finish_ = str_->size();
  if ((ssize_t) (finish_ - start_) + 1 < l_)  // Invalidate scanner if interval < 1 l-mer
    str_pos_ = finish_;
}

void MinimizerScanner::LoadSequence(string &seq, size_t start, size_t finish) {
  str_ = &seq;
  start_ = start;
  finish_ = finish;
  str_pos_ = start_;
  if (finish_ == SIZE_MAX)
    finish_ = str_->size();
  if ((ssize_t) (finish_ - start_) + 1 < l_)  // Invalidate scanner if interval < 1 l-mer
    str_pos_ = finish_;
  queue_.clear();
  queue_pos_ = 0;
  loaded_ch_ = 0;
  last_minimizer_ = ~0;
  last_ambig_ = 0;
}

uint64_t *MinimizerScanner::NextMinimizer(bool *ambig_flag) {
  if (str_pos_ >= finish_)  // Abort if we've exhausted string interval
    return nullptr;
  bool changed_minimizer = false;
  while (! changed_minimizer) {
    // Incorporate next character (and more if needed to fill l-mer)
    if (loaded_ch_ == l_)
      loaded_ch_--;
    while (loaded_ch_ < l_ && str_pos_ < finish_) {  // char loading loop
      loaded_ch_++;
      if (dna_) {
        lmer_ <<= BITS_PER_CHAR_DNA;
        last_ambig_ <<= BITS_PER_CHAR_DNA;
        switch ((*str_)[str_pos_++]) {
          case 'A' : case 'a' :                break;
          case 'C' : case 'c' : lmer_ |= 0x01; break;
          case 'G' : case 'g' : lmer_ |= 0x02; break;
          case 'T' : case 't' : lmer_ |= 0x03; break;
          default:  // ambig code, dump nt from l-mer and wipe queue
            queue_.clear();
            queue_pos_ = 0;
            lmer_ = 0;
            loaded_ch_ = 0;
            last_ambig_ |= 0x03;
        }
      }
      else {
        lmer_ <<= BITS_PER_CHAR_PRO;
        last_ambig_ <<= BITS_PER_CHAR_PRO;
        switch ((*str_)[str_pos_++]) {
          // Reduced alphabet uses 15-letter alphabet from AD Solis (2015),
          // in Proteins (doi:10.1002/prot.24936)
          // This means that the ambiguous codes B (N/D) and Z (E/Q) aren't
          // usable here because they span multiple groups.  Although J (I/L)
          // could be used because I & L are in the same group, we treat it as
          // we do B, Z, and X for consistency and to ease future changes to
          // the code base.
          case '*' :             // protein termination
          case 'U' : case 'u' :  // selenocysteine (rare)
          case 'O' : case 'o' :  // pyrrolysine (rare)
            // stop codons/rare amino acids, just use zero
            break;
          case 'A' : case 'a' :  // alanine
            lmer_ |= 0x01; break;
          case 'N' : case 'n' :  // asparagine
          case 'Q' : case 'q' :  // glutamine
          case 'S' : case 's' :  // serine
            lmer_ |= 0x02; break;
          case 'C' : case 'c' :  // cysteine
            lmer_ |= 0x03; break;
          case 'D' : case 'd' :  // aspartic acid
          case 'E' : case 'e' :  // glutamic acid
            lmer_ |= 0x04; break;
          case 'F' : case 'f' :  // phenylalanine
            lmer_ |= 0x05; break;
          case 'G' : case 'g' :  // glycine
            lmer_ |= 0x06; break;
          case 'H' : case 'h' :  // histidine
            lmer_ |= 0x07; break;
          case 'I' : case 'i' :  // isoleucine
          case 'L' : case 'l' :  // leucine
            lmer_ |= 0x08; break;
          case 'K' : case 'k' :  // lysine
            lmer_ |= 0x09; break;
          case 'P' : case 'p' :  // proline
            lmer_ |= 0x0a; break;
          case 'R' : case 'r' :  // arginine
            lmer_ |= 0x0b; break;
          case 'M' : case 'm' :  // methionine
          case 'V' : case 'v' :  // valine
            lmer_ |= 0x0c; break;
          case 'T' : case 't' :  // threonine
            lmer_ |= 0x0d; break;
          case 'W' : case 'w' :  // tryptophan
            lmer_ |= 0x0e; break;
          case 'Y' : case 'y' :  // tyrosine
            lmer_ |= 0x0f; break;
          default:  // ambig code, dump aa from l-mer and wipe queue
            queue_.clear();
            queue_pos_ = 0;
            lmer_ = 0;
            loaded_ch_ = 0;
            last_ambig_ |= 0x0f;
        }
      }
      lmer_ &= lmer_mask_;
      last_ambig_ &= lmer_mask_;
      // non-null arg means we need to do some special handling
      if (ambig_flag != nullptr) {
        *ambig_flag = !! last_ambig_;
        // If we haven't filled up first k-mer, don't return
        // Otherwise, if l-mer is incomplete, still return
        // If l-mer is complete, go through queue management
        // Ambig flag being non-null means we're expecting a return for each
        // k-mer, and if we don't have a full l-mer, there's no point in
        // queue management for minimizer calculation
        if (str_pos_ >= (size_t) k_ && loaded_ch_ < l_)
          return &last_minimizer_;
      }
    }  // end character loading loop
    if (loaded_ch_ < l_)  // Abort if we've exhausted string interval w/o
      return nullptr;     //   filling the l-mer
    uint64_t canonical_lmer = dna_ ? canonical_representation(lmer_, l_) : lmer_;
    if (spaced_seed_mask_)
      canonical_lmer &= spaced_seed_mask_;
    uint64_t candidate_lmer = canonical_lmer ^ toggle_mask_;
    if (k_ == l_) {  // Short-circuit queue work
      last_minimizer_ = candidate_lmer ^ toggle_mask_;
      return &last_minimizer_;
    }
    // Sliding window minimum calculation
    while (! queue_.empty() && queue_.back().candidate > candidate_lmer)
      queue_.pop_back();
    MinimizerData data = { candidate_lmer, queue_pos_ };
    if (queue_.empty() && queue_pos_ >= k_ - l_) {
      // Empty queue means front will change
      // Minimizer will change iff we've processed enough l-mers
      changed_minimizer = true;
    }
    queue_.push_back(data);
    // expire an l-mer not in the current window
    if (queue_.front().pos < queue_pos_ - k_ + l_) {
      queue_.erase(queue_.begin());
      // Change in front means minimizer changed
      changed_minimizer = true;
    }

    // Change from no minimizer (beginning of sequence/near ambig. char)
    if (queue_pos_ == k_ - l_)
      changed_minimizer = true;
    queue_pos_++;

    // (almost) Always return after each ch. when ambig flag is set
    if (ambig_flag != nullptr) {
      // Return only if we've read in at least one k-mer's worth of chars
      if (str_pos_ >= (size_t) k_) {
        // Force flag to true if we haven't read in a k-mer since last clearing
        // Must use <= due to increment of queue_pos_ just above here
        if (queue_pos_ <= k_ - l_)
          *ambig_flag = true;
        break;
      }
    }
  }  // end while ! changed_minimizer
  assert(! queue_.empty());
  last_minimizer_ = queue_.front().candidate ^ toggle_mask_;
  return &last_minimizer_;
}

// Adapted for 64-bit DNA use from public domain code at:
//   https://graphics.stanford.edu/~seander/bithacks.html#ReverseParallel
uint64_t MinimizerScanner::reverse_complement(uint64_t kmer, uint8_t n) {
  // Reverse bits (leaving bit pairs - nucleotides - intact)
  // swap consecutive pairs
  kmer = ((kmer & 0xCCCCCCCCCCCCCCCCUL) >> 2)
       | ((kmer & 0x3333333333333333UL) << 2);
  // swap consecutive nibbles
  kmer = ((kmer & 0xF0F0F0F0F0F0F0F0UL) >> 4)
       | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
  // swap consecutive bytes
  kmer = ((kmer & 0xFF00FF00FF00FF00UL) >> 8)
       | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
  // swap consecutive byte pairs
  kmer = ((kmer & 0xFFFF0000FFFF0000UL) >> 16)
       | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
  // swap halves of 64-bit word
  kmer = ( kmer >> 32 ) | ( kmer << 32);
  // Then complement
  return (~kmer) & ((1ull << (n * 2)) - 1);
}

uint64_t MinimizerScanner::canonical_representation(uint64_t kmer, uint8_t n) {
  uint64_t revcom = reverse_complement(kmer, n);
  return kmer < revcom ? kmer : revcom;
}

} // end namespace
