/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_KV_STORE_H_
#define KRAKEN2_KV_STORE_H_

#include "kraken2_headers.h"

#include <string>
#include <fstream>

namespace kraken2 {

typedef uint64_t hkey_t;
typedef uint32_t hvalue_t;

enum CellType {
  CompactHash32,
  CompactHash40,
  Unknown,
};

class KeyValueStore {
  public:
  virtual hvalue_t Get(hkey_t key) const = 0;
  virtual ~KeyValueStore() {}

  // static CellType GetKVCellType(std::string &filename);
};

CellType inline GetKVStoreCellType(std::string &filename) {

  size_t capacity;
  size_t size;
  size_t key_bits;
  size_t value_bits;

  std::ifstream ifs(filename);
  ifs.read((char *) &capacity, sizeof(capacity));
  ifs.read((char *) &size, sizeof(size));
  ifs.read((char *) &key_bits, sizeof(key_bits));
  ifs.read((char *)&value_bits, sizeof(value_bits));

  switch (key_bits + value_bits) {
  case 32:
    return CompactHash32;
  case 40:
    return CompactHash40;
  default:
    // Should not get here
    return Unknown;
  }
 }

uint64_t inline MurmurHash3(hkey_t key) {
  uint64_t k = (uint64_t) key;
  k ^= k >> 33;
  k *= 0xff51afd7ed558ccd;
  k ^= k >> 33;
  k *= 0xc4ceb9fe1a85ec53;
  k ^= k >> 33;
  return k;
}

}  // end namespace

#endif
