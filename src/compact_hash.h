/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_COMPACT_HASH_H_
#define KRAKEN2_COMPACT_HASH_H_

#include "kv_store.h"
#include "mmap_file.h"
#include "kraken2_headers.h"
#include "kraken2_data.h"
#include <sys/mman.h>
#include <cstdlib>
#include <fcntl.h>
#include <unistd.h>
#include <cerrno>

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

struct CompactHashCell40 {
  uint32_t a;
  uint8_t b;

  // inline hkey_t hashed_key(hkey_t value_bits) {
  //   size_t key_bits = 32 - value_bits - 8;
  //   size_t key_mask = (1 << key_bits) - 1;
  //   return a & key_mask;
  // }

  inline hkey_t hashed_key(size_t value_bits) {
    size_t key_bits = 32 - value_bits;
    return (a >> value_bits | b << key_bits);
  }

  inline hvalue_t value(size_t value_bits) {
    return (hvalue_t) (a & ((1 << value_bits) - 1));
  }


  void populate(hkey_t compacted_key, hvalue_t val, size_t key_bits, size_t value_bits) {
    if (key_bits + value_bits != sizeof(CompactHashCell40) * 8)
      errx(EX_SOFTWARE, "key len of %u and value len of %u don't sum to %d",
           (unsigned int) key_bits, (unsigned int) value_bits, (unsigned int)sizeof(CompactHashCell40));
    if (! key_bits || ! value_bits)
      errx(EX_SOFTWARE, "key len and value len must be nonzero");
    uint64_t max_value = (1llu << value_bits) - 1;
    if (max_value < val)
      errx(EX_SOFTWARE, "value len of %u too small for value of %llu",
           (unsigned int) value_bits, (unsigned long long int) val);

    // set the value
    size_t value_mask = (1 << value_bits) - 1;
    val = val & value_mask;

    // set the key
    size_t a_bits = 32 - value_bits;
    b = compacted_key >> a_bits;
    a = ((compacted_key & ((1 << a_bits) - 1)) << value_bits) | val;
  }
} __attribute__((packed));

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

template<typename Cell> class CompactHashTable : public KeyValueStore {
  public:
  CompactHashTable(size_t capacity, size_t key_bits, size_t value_bits);
  CompactHashTable(const std::string &filename, bool memory_mapping=false);
  CompactHashTable(const char *filename, bool memory_mapping=false);
  ~CompactHashTable();

  hvalue_t Get(hkey_t key) const;
  void GetBatch(const hkey_t *keys, hvalue_t *out, size_t n) const;
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
  Cell *table_;
  bool file_backed_;
  bool locks_initialized_;
  MMapFile backing_file_;
  omp_lock_t zone_locks_[LOCK_ZONES];

  CompactHashTable(const CompactHashTable &rhs);
  CompactHashTable& operator=(const CompactHashTable &rhs);

  void LoadTable(const char *filename, bool memory_mapping);
  uint64_t second_hash(uint64_t first_hash) const;
};

        template<typename Cell>
CompactHashTable<Cell>::CompactHashTable(size_t capacity, size_t key_bits, size_t value_bits)
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
  {
    size_t table_bytes = capacity_ * sizeof(*table_);
    if (posix_memalign((void **) &table_, (size_t) 1 << 21, table_bytes) != 0) {
      std::cerr << "Failed attempt to allocate " << table_bytes << "bytes;\n"
                << "you may not have enough free memory to build this database.\n"
                << "Perhaps increasing the k-mer length, or reducing memory usage from\n"
                << "other programs could help you build this database?" << std::endl;
      errx(EX_OSERR, "unable to allocate hash table memory");
    }
#ifdef MADV_HUGEPAGE
    madvise(table_, table_bytes, MADV_HUGEPAGE);
#endif
    memset(table_, 0, table_bytes);
  }
}

template<typename Cell>
CompactHashTable<Cell>::CompactHashTable(const string &filename, bool memory_mapping) {
  LoadTable(filename.c_str(), memory_mapping);
}

template<typename Cell>
CompactHashTable<Cell>::CompactHashTable(const char *filename, bool memory_mapping) {
  LoadTable(filename, memory_mapping);
}

template<typename Cell>
CompactHashTable<Cell>::~CompactHashTable() {
  if (! file_backed_)
    free(table_);
  if (locks_initialized_)
    for (size_t i = 0; i < LOCK_ZONES; i++)
      omp_destroy_lock(&zone_locks_[i]);
}

// Read exactly n bytes from fd into buf, looping over short reads (a single
// read() is capped well below multi-GB sizes), retrying on EINTR. errx on
// error or premature EOF.
static inline void read_fully(int fd, void *buf, size_t n, const char *what) {
  char *p = (char *) buf;
  while (n > 0) {
    ssize_t r = read(fd, p, n);
    if (r < 0) {
      if (errno == EINTR) continue;
      errx(EX_OSERR, "error reading %s", what);
    }
    if (r == 0)
      errx(EX_DATAERR, "unexpected end of file reading %s", what);
    p += r;
    n -= (size_t) r;
  }
}

// Number of concurrent reader threads for loading the table; tunable via
// K2_DB_READ_THREADS (default 8). More streams raise I/O queue depth, which is
// what saturates the device — reading is I/O-bound, not CPU-bound, so this can
// exceed core count.
static inline int db_read_threads() {
  const char *e = getenv("K2_DB_READ_THREADS");
  int n = e ? atoi(e) : 8;
  return n > 0 ? n : 1;
}

// Fill [buf, buf+n) from fd starting at file offset base_off using `threads`
// concurrent pread() streams over disjoint chunks. Each pread carries its own
// offset and is thread-safe. errx on error/EOF.
static inline void pread_parallel(int fd, void *buf, size_t n, off_t base_off,
                                  int threads, const char *what) {
  if (threads < 1) threads = 1;
  size_t chunk = (n + (size_t) threads - 1) / (size_t) threads;
  // Set to 1 by any thread that hits a read error. Plain assignment is safe here:
  // every writer stores the same value, and the read below runs after the parallel
  // for's implicit barrier (no concurrent access), so no atomics are needed -- and
  // it avoids OpenMP 3.1's `atomic write` for broader toolchain compatibility.
  int failed = 0;
  #pragma omp parallel for schedule(static) num_threads(threads)
  for (int t = 0; t < threads; t++) {
    size_t start = (size_t) t * chunk;
    if (start >= n) continue;
    size_t len = chunk < n - start ? chunk : n - start;
    char *p = (char *) buf + start;
    off_t off = base_off + (off_t) start;
    size_t got = 0;
    while (got < len) {
      ssize_t r = pread(fd, p + got, len - got, off + (off_t) got);
      if (r < 0) {
        if (errno == EINTR) continue;
        failed = 1;
        break;
      }
      if (r == 0) {  // premature EOF
        failed = 1;
        break;
      }
      got += (size_t) r;
    }
  }
  if (failed)
    errx(EX_OSERR, "error reading %s", what);
}

template<typename Cell>
void CompactHashTable<Cell>::LoadTable(const char *filename, bool memory_mapping) {
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
    table_ = (Cell *) ptr;
    if (backing_file_.filesize() - (ptr - backing_file_.fptr()) !=
        sizeof(*table_) * capacity_)
    {
      errx(EX_DATAERR, "Capacity mismatch in %s, aborting", filename);
    }
    file_backed_ = true;
  }
  else {
    int fd = open(filename, O_RDONLY);
    if (fd < 0)
      errx(EX_OSERR, "unable to open %s", filename);
    read_fully(fd, &capacity_, sizeof(capacity_), filename);
    read_fully(fd, &size_, sizeof(size_), filename);
    read_fully(fd, &key_bits_, sizeof(key_bits_), filename);
    read_fully(fd, &value_bits_, sizeof(value_bits_), filename);
    off_t table_off = lseek(fd, 0, SEEK_CUR);  // byte offset where the table begins
    if (table_off < 0)
      errx(EX_OSERR, "lseek failed on %s", filename);
    {
      size_t table_bytes = capacity_ * sizeof(*table_);
      // 2MB-aligned so MADV_HUGEPAGE (where supported) can back the multi-GB
      // table with huge pages, cutting TLB-miss page-walks on the random-access
      // lookup path (large win on aarch64/Graviton). Cell is POD; free() pairs
      // with this. On platforms without THP the madvise is simply compiled out.
      if (posix_memalign((void **) &table_, (size_t) 1 << 21, table_bytes) != 0) {
        std::cerr << "Failed attempt to allocate " << table_bytes << "bytes;\n"
                  << "you may not have enough free memory to load this database.\n"
                  << "If your computer has enough RAM, perhaps reducing memory usage from\n"
                  << "other programs could help you load this database?" << std::endl;
        errx(EX_OSERR, "unable to allocate hash table memory");
      }
#ifdef MADV_HUGEPAGE
      madvise(table_, table_bytes, MADV_HUGEPAGE);
#endif
      // Read the multi-GB table with several concurrent pread() streams. A single
      // sequential read keeps too few requests in flight (bounded by the small
      // default readahead window) to saturate a fast block device such as EBS, so
      // load is throttled well below device throughput. Splitting the table across
      // N threads, each pread()-ing a disjoint region, raises the effective I/O
      // queue depth and lets the device run at full speed. pread is thread-safe
      // (it does not use or move the shared file offset). This is a pure userspace
      // change: no privileges, no host/readahead tuning, identical in a container.
      // Without OpenMP the pragma is ignored and this degrades to a serial read.
      pread_parallel(fd, table_, table_bytes, table_off, db_read_threads(), filename);
    }
    close(fd);
    file_backed_ = false;
  }
}

template<typename Cell>
void CompactHashTable<Cell>::WriteTable(const char *filename) {
  ofstream ofs(filename, ofstream::binary);
  ofs.write((char *) &capacity_, sizeof(capacity_));
  ofs.write((char *) &size_, sizeof(size_));
  ofs.write((char *) &key_bits_, sizeof(key_bits_));
  ofs.write((char *) &value_bits_, sizeof(value_bits_));
  ofs.write((char *) table_, sizeof(*table_) * capacity_);
  ofs.close();
}

template<typename Cell>
hvalue_t CompactHashTable<Cell>::Get(hkey_t key) const {
  uint64_t hc = MurmurHash3(key);
  uint64_t compacted_key = hc >> (64 - key_bits_);
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

// Resolve a batch of keys, software-prefetching the initial probe cell of a key
// PREFETCH_DIST iterations ahead so multiple cache-missing loads are in flight at
// once (memory-level parallelism). One virtual dispatch per batch amortizes away.
template<typename Cell>
void CompactHashTable<Cell>::GetBatch(const hkey_t *keys, hvalue_t *out, size_t n) const {
  static const size_t PREFETCH_DIST = [] {
    const char *e = getenv("K2_PREFETCH_DIST");
    return e ? (size_t) atoi(e) : (size_t) 8;
  }();
  for (size_t i = 0; i < n; i++) {
    if (i + PREFETCH_DIST < n) {
      uint64_t hc = MurmurHash3(keys[i + PREFETCH_DIST]);
      __builtin_prefetch(&table_[hc % capacity_], 0, 1);
    }
    out[i] = Get(keys[i]);
  }
}

template<typename Cell>
bool CompactHashTable<Cell>::FindIndex(hkey_t key, size_t *idx) const {
  uint64_t hc = MurmurHash3(key);
  uint64_t compacted_key = hc >> (64 - key_bits_);
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

template<typename Cell>
bool CompactHashTable<Cell>::CompareAndSet
    (hkey_t key, hvalue_t new_value, hvalue_t *old_value)
{
  if (file_backed_)
    return false;
  if (new_value == 0)
    return false;
  uint64_t hc = MurmurHash3(key);
  hkey_t compacted_key = hc >> (64 - key_bits_);
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
    if (search_successful)
      break;
    if (step == 0)
      step = second_hash(hc);
    idx += step;
    idx %= capacity_;
    if (idx == first_idx)
      errx(EX_SOFTWARE, "compact hash table capacity exceeded");
  }
  return set_successful;
}

template<typename Cell>
bool CompactHashTable<Cell>::DirectCompareAndSet
    (size_t idx, hkey_t key, hvalue_t new_value, hvalue_t *old_value)
{
  uint64_t hc = MurmurHash3(key);
  hkey_t compacted_key = hc >> (64 - key_bits_);
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
template<typename Cell>
inline uint64_t CompactHashTable<Cell>::second_hash(uint64_t first_hash) const {
#ifdef LINEAR_PROBING
  return 1;
#else  // Double hashing
  return (first_hash >> 8) | 1;
#endif
}

template<typename Cell>
taxon_counts_t CompactHashTable<Cell>::GetValueCounts() const {
  taxon_counts_t value_counts;
  int thread_ct = omp_get_max_threads();
  taxon_counts_t *thread_value_counts = new taxon_counts_t[thread_ct];
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
  delete[] thread_value_counts;
  return value_counts;
}


}  // end namespace

#endif
