/*
 * Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_KV_STORE_H_
#define KRAKEN2_KV_STORE_H_

#include "kraken2_headers.h"

namespace kraken2 {

typedef uint64_t hkey_t;
typedef uint32_t hvalue_t;

class KeyValueStore {
  public:
  virtual hvalue_t Get(hkey_t key) const = 0;
  virtual ~KeyValueStore() { }
};

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
