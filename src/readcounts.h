/*
 * Copyright 2017-2018, Florian Breitwieser
 *
 * This file was originally developed for the KrakenUniq taxonomic classification system.
 */

// Note: This class provides counter information for read and k-mer data
// as well as allowing counting/estimation of distinct k-mers via a container
// type that is passed in.

#ifndef KRAKEN2_READCOUNTS_H
#define KRAKEN2_READCOUNTS_H

#include "hyperloglogplus.h"
#include <unordered_set>

namespace kraken2 {

  template <typename CONTAINER>
  class ReadCounts {

  public:
    uint64_t readCount() const { return n_reads; }
    void incrementReadCount() { ++n_reads; }
    uint64_t kmerCount() const { return n_kmers; }
    uint64_t distinctKmerCount() const; // to be implemented for each CONTAINER

    ReadCounts() : n_reads(0), n_kmers(0) {
    }

    ReadCounts(uint64_t _n_reads, uint64_t _n_kmers) : n_reads(_n_reads), n_kmers(_n_kmers) {
    }

    ReadCounts(const ReadCounts& other) : n_reads(other.n_reads), n_kmers(other.n_kmers), kmers(other.kmers) {
    }

    ReadCounts(ReadCounts&& other) : n_reads(other.n_reads), n_kmers(other.n_kmers), kmers(std::move(other.kmers)) {
    }

    ReadCounts& operator=(const ReadCounts& other) {
      n_reads = other.n_reads;
      n_kmers = other.n_kmers;
      kmers = other.kmers;
      return *this;
    }


    ReadCounts& operator=(ReadCounts&& other) {
      n_reads = other.n_reads;
      n_kmers = other.n_kmers;
      kmers = std::move(other.kmers);
      return *this;
    }

    void add_kmer(uint64_t kmer) {
      ++n_kmers;
      kmers.insert(kmer);
    }

    ReadCounts& operator+=(const ReadCounts& other) {
      n_reads += other.n_reads;
      n_kmers += other.n_kmers;
      kmers += other.kmers;
      return *this;
    }

    ReadCounts& operator+=(ReadCounts&& other) {
      n_reads += other.n_reads;
      n_kmers += other.n_kmers;
      kmers += std::move(other.kmers);
      return *this;
    }

    bool operator<(const ReadCounts& other) {
      if (n_reads < other.n_reads) {
        return true;
      }
      if (n_reads == other.n_reads && n_kmers < other.n_kmers) {
        return true;
      }
      return false;
    }

  private:
    uint64_t n_reads;
    uint64_t n_kmers;
    CONTAINER kmers; // distinct k-mer count per taxon
  };

  // Overload operator += for set, so that it can be used for merging
  template <typename T>
  unordered_set<T>& operator+=(unordered_set<T>& left, const unordered_set<T>& right) {
    left.insert(right.begin(), right.end());
    return left;
  }

  template <typename T>
  set<T>& operator+=(set<T>& left, const set<T>& right) {
    left.insert(right.begin(), right.end());
    return left;
  }

  template<typename T>
  uint64_t ReadCounts< T >::distinctKmerCount() const {
    // DEW: changed cardinality() to size() to enable use of STL containers
    return kmers.size();
  }

}
#endif
