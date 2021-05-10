/*
 * Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include "compact_hash.h"

using std::string;
using std::ofstream;
using std::unordered_map;

namespace kraken2 {

CompactHashTable::CompactHashTable(size_t capacity, size_t key_bits, size_t value_bits)
    : capacity_(capacity), size_(0), key_bits_(key_bits), value_bits_(value_bits),
      file_backed_(false), locks_initialized_(true)
{
  if (key_bits + value_bits != sizeof(*table_) * 8)
    errx(EX_SOFTWARE, "sum of key bits and value bits must equal %u",
        (unsigned int) (sizeof(*table_) * 8));
  if (key_bits == 0)
    errx(EX_SOFTWARE, "key bits cannot be zero");
  if (value_bits == 0)
    errx(EX_SOFTWARE, "value bits cannot be zero");
  for (size_t i = 0; i < LOCK_ZONES; i++)
    omp_init_lock(&zone_locks_[i]);
  try {
    table_ = new CompactHashCell[capacity_];
  } catch (std::bad_alloc ex) {
    std::cerr << "Failed attempt to allocate " << (sizeof(*table_) * capacity_) << "bytes;\n"
              << "you may not have enough free memory to build this database.\n"
              << "Perhaps increasing the k-mer length, or reducing memory usage from\n"
              << "other programs could help you build this database?" << std::endl;
    errx(EX_OSERR, "unable to allocate hash table memory");
  }
  memset(table_, 0, capacity_ * sizeof(*table_));
}

CompactHashTable::CompactHashTable(const string &filename, bool memory_mapping) {
  LoadTable(filename.c_str(), memory_mapping);
}

CompactHashTable::CompactHashTable(const char *filename, bool memory_mapping) {
  LoadTable(filename, memory_mapping);
}

CompactHashTable::~CompactHashTable() {
  if (! file_backed_)
    delete[] table_;
  if (locks_initialized_)
    for (size_t i = 0; i < LOCK_ZONES; i++)
      omp_destroy_lock(&zone_locks_[i]);
}

void CompactHashTable::LoadTable(const char *filename, bool memory_mapping) {
  locks_initialized_ = false;
  if (memory_mapping) {
    backing_file_.OpenFile(filename);
    char *ptr = backing_file_.fptr();
    memcpy((char *) &capacity_, ptr, sizeof(capacity_));
    ptr += sizeof(capacity_);
    memcpy((char *) &size_, ptr, sizeof(size_));
    ptr += sizeof(size_);
    memcpy((char *) &key_bits_, ptr, sizeof(key_bits_));
    ptr += sizeof(key_bits_);
    memcpy((char *) &value_bits_, ptr, sizeof(value_bits_));
    ptr += sizeof(value_bits_);
    table_ = (CompactHashCell *) ptr;
    if (backing_file_.filesize() - (ptr - backing_file_.fptr()) !=
        sizeof(*table_) * capacity_)
    {
      errx(EX_DATAERR, "Capacity mismatch in %s, aborting", filename);
    }
    file_backed_ = true;
  }
  else {
    std::ifstream ifs(filename);
    ifs.read((char *) &capacity_, sizeof(capacity_));
    ifs.read((char *) &size_, sizeof(size_));
    ifs.read((char *) &key_bits_, sizeof(key_bits_));
    ifs.read((char *) &value_bits_, sizeof(value_bits_));
    try {
      table_ = new CompactHashCell[capacity_];
    } catch (std::bad_alloc &ex) {
      std::cerr << "Failed attempt to allocate " << (sizeof(*table_) * capacity_) << "bytes;\n"
                << "you may not have enough free memory to load this database.\n"
                << "If your computer has enough RAM, perhaps reducing memory usage from\n"
                << "other programs could help you load this database?" << std::endl;
      errx(EX_OSERR, "unable to allocate hash table memory");
    }
    ifs.read((char *) table_, capacity_ * sizeof(*table_));
    if (! ifs)
      errx(EX_OSERR, "Error reading in hash table");
    file_backed_ = false;
  }
}

void CompactHashTable::WriteTable(const char *filename) {
  ofstream ofs(filename, ofstream::binary);
  ofs.write((char *) &capacity_, sizeof(capacity_));
  ofs.write((char *) &size_, sizeof(size_));
  ofs.write((char *) &key_bits_, sizeof(key_bits_));
  ofs.write((char *) &value_bits_, sizeof(value_bits_));
  ofs.write((char *) table_, sizeof(*table_) * capacity_);
  ofs.close();
}

hvalue_t CompactHashTable::Get(hkey_t key) const {
  uint64_t hc = MurmurHash3(key);
  uint64_t compacted_key = hc >> (32 + value_bits_);
  size_t idx = hc % capacity_;
  size_t first_idx = idx;
  size_t step = 0;
  while (true) {
    if (! table_[idx].value(value_bits_))  // value of 0 means data is 0, saves work
      break;  // search over, empty cell encountered in probe
    if (table_[idx].hashed_key(value_bits_) == compacted_key)
      return table_[idx].value(value_bits_);
    if (step == 0)
      step = second_hash(hc);
    idx += step;
    idx %= capacity_;
    if (idx == first_idx)
      break;  // search over, we've exhausted the table
  }
  return 0;
}

bool CompactHashTable::FindIndex(hkey_t key, size_t *idx) const {
  uint64_t hc = MurmurHash3(key);
  uint64_t compacted_key = hc >> (32 + value_bits_);
  *idx = hc % capacity_;
  size_t first_idx = *idx;
  size_t step = 0;
  while (true) {
    if (! table_[*idx].value(value_bits_))  // value of 0 means data is 0, saves work
      return false;  // search over, empty cell encountered in probe
    if (table_[*idx].hashed_key(value_bits_) == compacted_key)
      return true;
    if (step == 0)
      step = second_hash(hc);
    *idx += step;
    *idx %= capacity_;
    if (*idx == first_idx)
      break;  // search over, we've exhausted the table
  }
  return false;
}

bool CompactHashTable::CompareAndSet
    (hkey_t key, hvalue_t new_value, hvalue_t *old_value)
{
  if (file_backed_)
    return false;
  if (new_value == 0)
    return false;
  uint64_t hc = MurmurHash3(key);
  hkey_t compacted_key = hc >> (32 + value_bits_);
  size_t idx, first_idx;
  bool set_successful = false;
  bool search_successful = false;
  idx = first_idx = hc % capacity_;
  size_t step = 0;
  while (! search_successful) {
    size_t zone = idx % LOCK_ZONES;
    omp_set_lock(&zone_locks_[zone]);
    if (! table_[idx].value(value_bits_)
        || table_[idx].hashed_key(value_bits_) == compacted_key)
    {
      search_successful = true;
      #pragma omp flush
      if (*old_value == table_[idx].value(value_bits_)) {
        table_[idx].populate(compacted_key, new_value, key_bits_, value_bits_);
        if (! *old_value) {
          #pragma omp atomic
          size_++;
        }
        set_successful = true;
      }
      else {
        *old_value = table_[idx].value(value_bits_);
      }
    }
    omp_unset_lock(&zone_locks_[zone]);
    if (step == 0)
      step = second_hash(hc);
    idx += step;
    idx %= capacity_;
    if (idx == first_idx)
      errx(EX_SOFTWARE, "compact hash table capacity exceeded");
  }
  return set_successful;
}

bool CompactHashTable::DirectCompareAndSet
    (size_t idx, hkey_t key, hvalue_t new_value, hvalue_t *old_value)
{
  uint64_t hc = MurmurHash3(key);
  hkey_t compacted_key = hc >> (32 + value_bits_);
  bool set_successful = false;
  size_t zone = idx % LOCK_ZONES;
  omp_set_lock(&zone_locks_[zone]);
  #pragma omp flush
  if (*old_value == table_[idx].value(value_bits_)) {
    table_[idx].populate(compacted_key, new_value, key_bits_, value_bits_);
    if (! *old_value) {
      #pragma omp atomic
      size_++;
    }
    set_successful = true;
  }
  else {
    *old_value = table_[idx].value(value_bits_);
  }
  omp_unset_lock(&zone_locks_[zone]);
  return set_successful;
}

// Linear probing may be ok for accuracy, as long as occupancy is < 95%
// Linear probing leads to more clustering, longer probing paths, and
//   higher probability of a false answer
// Double hashing can have shorter probing paths, but less cache efficiency
inline uint64_t CompactHashTable::second_hash(uint64_t first_hash) const {
#ifdef LINEAR_PROBING
  return 1;
#else  // Double hashing
  return (first_hash >> 8) | 1;
#endif
}

taxon_counts_t CompactHashTable::GetValueCounts() const {
  taxon_counts_t value_counts;
  int thread_ct = omp_get_max_threads();
  taxon_counts_t thread_value_counts[thread_ct];
  #pragma omp parallel for
  for (size_t i = 0; i < capacity_; i++) {
    auto val = table_[i].value(value_bits_);
    if (val)
      thread_value_counts[omp_get_thread_num()][val]++;
  }
  for (auto i = 0; i < thread_ct; i++) {
    for (auto &kv_pair : thread_value_counts[i])
      value_counts[kv_pair.first] += kv_pair.second;
  }
  return value_counts;
}

}  // end namespace
