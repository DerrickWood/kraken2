/*
 * Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_COMPACT_HASH_H_
#define KRAKEN2_COMPACT_HASH_H_

#include "kv_store.h"
#include "mmap_file.h"
#include "kraken2_headers.h"
#include "kraken2_data.h"

namespace kraken2 {

struct CompactHashCell {
  inline hkey_t hashed_key(size_t value_bits) {
    return (hkey_t) (data >> value_bits);
  }

  inline hvalue_t value(size_t value_bits) {
    return (hvalue_t) (data & ((1 << value_bits) - 1));
  }

  void populate(hkey_t compacted_key, hvalue_t val, size_t key_bits, size_t value_bits) {
    if (key_bits + value_bits != 32)
      errx(EX_SOFTWARE, "key len of %u and value len of %u don't sum to 32",
           (unsigned int) key_bits, (unsigned int) value_bits);
    if (! key_bits || ! value_bits)
      errx(EX_SOFTWARE, "key len and value len must be nonzero");
    uint64_t max_value = (1llu << value_bits) - 1;
    if (max_value < val)
      errx(EX_SOFTWARE, "value len of %u too small for value of %llu",
           (unsigned int) value_bits, (unsigned long long int) val);
    data = compacted_key << value_bits;
    data |= val;
  }

  // value in the low bits
  // hash of the key in the high bits
  uint32_t data;
};

/**
 Implementation of a compact, fixed-size probabilistic hash table.

 Operations supported:

 - search
 - compare-and-set

 Does NOT support:

 - values of 0

 Keys must be 64-bit unsigned integers, < 2 ** 63.
 Values must be 32-bit unsigned integers, > 0.

 Cells are 32-bit unsigned integers, with truncated hashed keys stored in
 the most significant bits, and truncated values stored in the least
 significant bits.  Compaction level is set when the table is created
 by the WriteCompactedTable() method of the HashTable class.
 **/

class CompactHashTable : public KeyValueStore {
  public:
  CompactHashTable(size_t capacity, size_t key_bits, size_t value_bits);
  CompactHashTable(const std::string &filename, bool memory_mapping=false);
  CompactHashTable(const char *filename, bool memory_mapping=false);
  ~CompactHashTable();

  hvalue_t Get(hkey_t key) const;
  bool FindIndex(hkey_t key, size_t *idx) const;

  // How CompareAndSet works:
  // if *old_value == CHT[key]
  //   CHT[key] = new_value
  //   return true
  // else
  //   *old_value = CHT[key]
  //   return false
  bool CompareAndSet(hkey_t key, hvalue_t new_value, hvalue_t *old_value);
  bool DirectCompareAndSet(size_t idx, hkey_t key, hvalue_t new_value, hvalue_t *old_value);
  void WriteTable(const char *filename);

  taxon_counts_t GetValueCounts() const;

  size_t capacity() const { return capacity_; }
  size_t size() const { return size_; }
  size_t key_bits() const { return key_bits_; }
  size_t value_bits() const { return value_bits_; }
  double occupancy() const { return size_ * 1.0 / capacity_; }

  private:
  static const size_t LOCK_ZONES = 256;
  size_t capacity_;
  size_t size_;
  //size_t index_mask_;
  size_t key_bits_;
  size_t value_bits_;
  CompactHashCell *table_;
  bool file_backed_;
  bool locks_initialized_;
  MMapFile backing_file_;
  omp_lock_t zone_locks_[LOCK_ZONES];

  CompactHashTable(const CompactHashTable &rhs);
  CompactHashTable& operator=(const CompactHashTable &rhs);

  void LoadTable(const char *filename, bool memory_mapping);
  uint64_t second_hash(uint64_t first_hash) const;
};

}  // end namespace

#endif
