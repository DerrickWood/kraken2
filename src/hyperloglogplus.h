/*
 * Copyright 2017-2018, Florian Breitwieser
 *
 * This file was originally developed for the KrakenUniq taxonomic classification system.
 */

/*
 * Implementation of 64-bit HyperLogLog algorithm by Flajolet et al.,
 *   with sparse mode for increased precision with low cardinalities (Stefan Heule et al.), and
 *   an improved estimator that does not rely on empirical bias correction data (Otmar Ertl)
 */

#ifndef HYPERLOGLOGPLUS_H_
#define HYPERLOGLOGPLUS_H_

#include "kraken2_headers.h"
using namespace std;

//#define HLL_DEBUG
//#define HLL_DEBUG2

#ifdef HLL_DEBUG
#define D(x) x
#else
#define D(x)
#endif

// 64-bit Mixer/finalizer from MurMurHash3 https://github.com/aappleby/smhasher
uint64_t murmurhash3_finalizer (uint64_t key);


// Heule et al. encode the sparse list with variable length encoding
//   see section 5.3.2. This implementation just uses a sorted vector or unordered_set.
typedef unordered_set<uint32_t> SparseListType;
// Other possible SparseList types:
// // typedef vector<uint32_t> SparseListType;
// The sorted vector implementation is pretty inefficient currently, as the vector
// is always kept sorted.
// //typedef unordered_set<uint32_t, NoHash<uint32_t> > SparseListType;
// No real performance gain from using the hash value directly - and probably problems
// with the bucket assignment, since unordered_set expects size_t hashes

/**
 * HyperLogLogPlusMinus class for counting the number of unique 64-bit values in stream
 * Note that only HASH=uint64_t is implemented.
 */
template<typename HASH>
class HyperLogLogPlusMinus {

private:
  uint8_t p;      // precision, set in constructor
  size_t m = 1 << p;  // number of registers
  vector<uint8_t> M;    // registers, size m
  uint64_t n_observed = 0;

  bool sparse;          // sparse representation of the data?
  SparseListType sparseList;
  HASH (*bit_mixer) (uint64_t);

  // sparse versions of p and m
  static const uint8_t  pPrime = 25; // precision when using a sparse representation
                                     // which encodes the rank + index in 32 bits:
                                     //  25 bits for index +
                                     //   6 bits for rank +
                                     //   1 flag bit indicating if bits p..pPrime are 0
  static const uint32_t mPrime = 1 << pPrime; // 2^pPrime

public:
  bool use_n_observed = true; // return min(estimate, n_observed) instead of estimate

  // Construct HLL with precision bits
  HyperLogLogPlusMinus(uint8_t precision=12, bool sparse=true, HASH (*bit_mixer) (uint64_t) = murmurhash3_finalizer);
  HyperLogLogPlusMinus(const HyperLogLogPlusMinus<HASH>& other);
  HyperLogLogPlusMinus(HyperLogLogPlusMinus<HASH>&& other);
  HyperLogLogPlusMinus<HASH>& operator= (HyperLogLogPlusMinus<HASH>&& other);
  HyperLogLogPlusMinus<HASH>& operator= (const HyperLogLogPlusMinus<HASH>& other);
  ~HyperLogLogPlusMinus() {};
  void reset(); // Note: sets sparse=true

  // Add items or other HLL to this sketch
  void insert(uint64_t item);
  void insert(const vector<uint64_t>& items);

  // Merge another sketch into this one
  // TODO: assumes equal bit_mixers! but does not check that
  void merge(HyperLogLogPlusMinus<HASH>&& other);
  void merge(const HyperLogLogPlusMinus<HASH>& other);
  HyperLogLogPlusMinus<HASH>& operator+=(HyperLogLogPlusMinus<HASH>&& other);
  HyperLogLogPlusMinus<HASH>& operator+=(const HyperLogLogPlusMinus<HASH>& other);

  // Calculate cardinality estimates
  uint64_t cardinality() const; // returns ertlCardinality()
  uint64_t size() const; // returns ertlCardinality()
  // HLL++ estimator of Heule et al., 2015. Uses empirical bias correction factors
  uint64_t heuleCardinality(bool correct_bias = true) const;
  // Improved estimator of Ertl, 2017. Does not rely on empirical data
  uint64_t ertlCardinality() const;
  // Flajolet's cardinality estimator without bias correction
  uint64_t flajoletCardinality(bool use_sparse_precision = true) const;

  uint64_t nObserved() const;

private:
  void switchToNormalRepresentation();
  void addToRegisters(const SparseListType &sparseList);

};

#endif /* HYPERLOGLOGPLUS_H_ */
