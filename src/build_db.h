#ifndef KRAKEN2_BUILD_DB_H_
#define KRAKEN2_BUILD_DB_H_

#include "compact_hash.h"
#include "taxonomy.h"
#include "kv_store.h"
#include "seqreader.h"
#include "mmscanner.h"

using namespace kraken2;

struct Options {
  string ID_to_taxon_map_filename;
  string ncbi_taxonomy_directory;
  string hashtable_filename;
  string options_filename;
  string taxonomy_filename;
  size_t block_size;
  size_t subblock_size;
  size_t cht_cell_size;
  size_t requested_bits_for_taxid;
  int num_threads;
  bool input_is_protein;
  ssize_t k, l;
  size_t capacity;
  size_t maximum_capacity;
  uint64_t spaced_seed_mask;
  uint64_t toggle_mask;
  uint64_t min_clear_hash_value;
  bool deterministic_build;
  size_t cht_cell_bits;
};

vector<string> ExtractNCBISequenceIDs(const string &header);
template<typename Cell>
void ProcessSequenceFast(const string &seq, taxid_t taxid,
    CompactHashTable<Cell> &hash, const Taxonomy &tax, MinimizerScanner &scanner,
    uint64_t min_clear_hash_value);
template<typename Cell>
void ProcessSequence(const Options &opts, const string &seq, taxid_t taxid,
    CompactHashTable<Cell> &hash, const Taxonomy &tax);
template<typename Cell>
void ProcessSequencesFast(const Options &opts,
    const map<string, taxid_t> &ID_to_taxon_map,
    CompactHashTable<Cell> &kraken_index, const Taxonomy &taxonomy);
template<typename Cell>
void ProcessSequences(const Options &opts,
    const map<string, taxid_t> &ID_to_taxon_map,
    CompactHashTable<Cell> &kraken_index, const Taxonomy &taxonomy);
template<typename Cell>
void SetMinimizerLCA(CompactHashTable<Cell> &hash, uint64_t minimizer, taxid_t taxid,
    const Taxonomy &tax);
template <typename Cell>
void build(const Taxonomy &tax, map<string, taxid_t> &ID_to_taxon_map, const Options &opts, size_t actual_capacity, size_t bits_for_taxid);


template <typename Cell>
void build(const Taxonomy &taxonomy, map<string, taxid_t> &ID_to_taxon_map,
           const Options &opts, size_t actual_capacity, size_t bits_for_taxid) {
  size_t total_bits = sizeof(Cell) * 8;
  CompactHashTable<Cell> kraken_index(actual_capacity, total_bits - bits_for_taxid,
      bits_for_taxid);
  std::cerr << "CHT created with " << bits_for_taxid << " bits reserved for taxid." << std::endl;

  if (opts.deterministic_build)
    ProcessSequences(opts, ID_to_taxon_map, kraken_index, taxonomy);
  else
    ProcessSequencesFast(opts, ID_to_taxon_map, kraken_index, taxonomy);

  std::cerr << "Writing data to disk... " << std::flush;
  kraken_index.WriteTable(opts.hashtable_filename.c_str());
}

// A quick but nondeterministic build
template<typename Cell>
void ProcessSequencesFast(const Options &opts,
                          const map<string, taxid_t> &ID_to_taxon_map,
                          CompactHashTable<Cell> &kraken_index,
                          const Taxonomy &taxonomy) {
  size_t processed_seq_ct = 0;
  size_t processed_ch_ct = 0;

#pragma omp parallel
  {
    Sequence sequence;
    MinimizerScanner scanner(opts.k, opts.l, opts.spaced_seed_mask,
                             !opts.input_is_protein, opts.toggle_mask);

    BatchSequenceReader reader;

    while (true) {
      // Declaration of "ok" and break need to be done outside of critical
      // section to conform with OpenMP spec.
      bool ok;
      BatchSequenceReader reader_clone(reader);
#pragma omp critical(reader)
      ok = reader_clone.LoadBlock(opts.block_size);
      if (!ok)
        break;
      while (reader_clone.NextSequence(sequence)) {
        auto all_sequence_ids = ExtractNCBISequenceIDs(sequence.header);
        taxid_t taxid = 0;
        for (auto &seqid : all_sequence_ids) {
          if (ID_to_taxon_map.count(seqid) == 0 ||
              ID_to_taxon_map.at(seqid) == 0)
            continue;
          auto ext_taxid = ID_to_taxon_map.at(seqid);
          taxid = taxonomy.LowestCommonAncestor(
              taxid, taxonomy.GetInternalID(ext_taxid));
        }
        if (taxid) {
          // Add terminator for protein sequences if not already there
          if (opts.input_is_protein && sequence.seq.back() != '*')
            sequence.seq.push_back('*');
          ProcessSequenceFast(sequence.seq, taxid, kraken_index, taxonomy,
                              scanner, opts.min_clear_hash_value);
#pragma omp atomic
          processed_seq_ct++;
#pragma omp atomic
          processed_ch_ct += sequence.seq.size();
        }
      }
      if (isatty(fileno(stderr))) {
#pragma omp critical(status_update)
        std::cerr << "\rProcessed " << processed_seq_ct << " sequences ("
                  << processed_ch_ct << " "
                  << (opts.input_is_protein ? "aa" : "bp") << ")...";
      }
    }
  }
  if (isatty(fileno(stderr)))
    std::cerr << "\r";
  std::cerr << "Completed processing of " << processed_seq_ct << " sequences, "
            << processed_ch_ct << " " << (opts.input_is_protein ? "aa" : "bp")
            << std::endl;
}

// Slightly slower but deterministic when multithreaded
template<typename Cell>
void ProcessSequences(const Options &opts,
                      const map<string, taxid_t> &ID_to_taxon_map,
                      CompactHashTable<Cell> &kraken_index,
                      const Taxonomy &taxonomy) {
  size_t processed_seq_ct = 0;
  size_t processed_ch_ct = 0;

  Sequence *sequence;
  BatchSequenceReader reader;

  while (reader.LoadBlock(opts.block_size)) {
    while ((sequence = reader.NextSequence()) != NULL) {
      auto all_sequence_ids = ExtractNCBISequenceIDs(sequence->header);
      taxid_t taxid = 0;
      int ext_taxid;
      for (auto &seqid : all_sequence_ids) {
        if (ID_to_taxon_map.count(seqid) == 0 ||
            ID_to_taxon_map.at(seqid) == 0) {
          continue;
        }
        ext_taxid = ID_to_taxon_map.at(seqid);
        taxid = taxonomy.LowestCommonAncestor(
            taxid, taxonomy.GetInternalID(ext_taxid));
      }
      if (taxid) {
        // Add terminator for protein sequences if not already there
        if (opts.input_is_protein && sequence->seq.back() != '*')
          sequence->seq.push_back('*');
        ProcessSequence(opts, sequence->seq, taxid, kraken_index, taxonomy);
        processed_seq_ct++;
        processed_ch_ct += sequence->seq.size();
      }
    }
    if (isatty(fileno(stderr))) {
      std::cerr << "\rProcessed " << processed_seq_ct << " sequences ("
                << processed_ch_ct << " "
                << (opts.input_is_protein ? "aa" : "bp") << ")...";
    }
  }
  if (isatty(fileno(stderr)))
    std::cerr << "\r";
  std::cerr << "Completed processing of " << processed_seq_ct << " sequences, "
            << processed_ch_ct << " " << (opts.input_is_protein ? "aa" : "bp")
            << std::endl;
}

template<typename Cell>
void SetMinimizerLCA(CompactHashTable<Cell> &hash, uint64_t minimizer, taxid_t taxid,
    const Taxonomy &tax)
{
  hvalue_t old_value = 0;
  hvalue_t new_value = taxid;
  while (! hash.CompareAndSet(minimizer, new_value, &old_value))
    new_value = tax.LowestCommonAncestor(old_value, taxid);
}

template<typename Cell>
void ProcessSequenceFast(const string &seq, taxid_t taxid,
    CompactHashTable<Cell> &hash, const Taxonomy &tax, MinimizerScanner &scanner,
    uint64_t min_clear_hash_value)
{
  scanner.LoadSequence(seq);
  uint64_t *minimizer_ptr;
  while ((minimizer_ptr = scanner.NextMinimizer())) {
    if (scanner.is_ambiguous())
      continue;
    if (min_clear_hash_value && MurmurHash3(*minimizer_ptr) < min_clear_hash_value)
      continue;
    hvalue_t existing_taxid = 0;
    hvalue_t new_taxid = taxid;
    while (! hash.CompareAndSet(*minimizer_ptr, new_taxid, &existing_taxid)) {
      new_taxid = tax.LowestCommonAncestor(new_taxid, existing_taxid);
    }
  }
}

template<typename Cell>
void ProcessSequence(const Options &opts, const string &seq, taxid_t taxid,
    CompactHashTable<Cell> &hash, const Taxonomy &tax)
{
  const int set_ct = 256;
  omp_lock_t locks[set_ct];
  for (int i = 0; i < set_ct; i++)
    omp_init_lock(&locks[i]);
  // for each block in the sequence
  for (size_t j = 0; j < seq.size(); j += opts.block_size) {
    // block: fixed length subsequence of DNA/protein, handled in series,
    //   consecutive blocks overlap by k-1 characters
    size_t block_start = j;
    size_t block_finish = j + opts.block_size + opts.k - 1;
    if (block_finish > seq.size())
      block_finish = seq.size();
    std::set<uint64_t> minimizer_sets[set_ct];

    // for each subblock in the block, gather minimizers in parallel
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = block_start; i < block_finish; i += opts.subblock_size) {
      // subblock: fixed length subsequence of block, handled in parallel,
      //   consecutive subblocks overlap by k-1 characters
      MinimizerScanner scanner(opts.k, opts.l, opts.spaced_seed_mask,
                               ! opts.input_is_protein, opts.toggle_mask);
      size_t subblock_finish = i + opts.subblock_size + opts.k - 1;
      if (subblock_finish > block_finish)
        subblock_finish = block_finish;
      scanner.LoadSequence(seq, i, subblock_finish);
      uint64_t *minimizer_ptr;
      while ((minimizer_ptr = scanner.NextMinimizer())) {
        if (scanner.is_ambiguous())
          continue;
        auto hc = MurmurHash3(*minimizer_ptr);
        // Hash-based subsampling
        if (opts.min_clear_hash_value && hc < opts.min_clear_hash_value)
          continue;
        auto zone = hc % set_ct;
        omp_set_lock(&locks[zone]);
        minimizer_sets[zone].insert(*minimizer_ptr);
        omp_unset_lock(&locks[zone]);
      }
    }  // end subblock for loop

    // combine sets into sorted list
    std::vector<uint64_t> minimizer_lists[set_ct];
    size_t minimizer_list_prefix_sizes[set_ct+1] = {0};
    #pragma omp parallel for
    for (int i = 0; i < set_ct; i++) {
      minimizer_lists[i].reserve(minimizer_sets[i].size());
      for (auto &m : minimizer_sets[i])
        minimizer_lists[i].push_back(m);
      sort(minimizer_lists[i].begin(), minimizer_lists[i].end());
      minimizer_list_prefix_sizes[i+1] = minimizer_lists[i].size();
    }
    for (int i = 2; i <= set_ct; i++)
      minimizer_list_prefix_sizes[i] += minimizer_list_prefix_sizes[i-1];
    std::vector<uint64_t> minimizer_list(minimizer_list_prefix_sizes[set_ct]);
    #pragma omp parallel for
    for (int i = 0; i < set_ct; i++) {
      std::copy(minimizer_lists[i].begin(), minimizer_lists[i].end(),
          minimizer_list.begin() + minimizer_list_prefix_sizes[i]);
    }

    size_t mm_ct = minimizer_list.size();
    std::vector<size_t> index_list(mm_ct);
    std::vector<bool> insertion_list(mm_ct);
    for (size_t i = 0; i < mm_ct; i++)
      insertion_list[i] = true;
    // Loop to enforce deterministic order
    while (! minimizer_list.empty()) {
      // Gather insertion point information for all remaining minimizers
      #pragma omp parallel for
      for (size_t i = 0; i < mm_ct; i++) {
        // once we've determined that a minimizer won't be an insertion,
        // don't bother calling FindIndex() again
        if (insertion_list[i]) {
          size_t idx;
          insertion_list[i] = ! hash.FindIndex(minimizer_list[i], &idx);
          index_list[i] = idx;
        }
      }

      // Determine safe prefix of sorted set to insert in parallel
      std::set<uint64_t> novel_insertion_points;
      size_t safe_ct = 0;
      for (safe_ct = 0; safe_ct < mm_ct; safe_ct++) {
        if (insertion_list[safe_ct]) {
          if (novel_insertion_points.count(index_list[safe_ct]) > 0)
            break;
          novel_insertion_points.insert(index_list[safe_ct]);
        }
      }

      // Adjust CHT values for all keys in the safe zone in parallel
      #pragma omp parallel for
      for (size_t i = 0; i < safe_ct; i++) {
        SetMinimizerLCA(hash, minimizer_list[i], taxid, tax);
      }

      // Remove safe prefix and re-iterate to process remainder
      minimizer_list.erase(minimizer_list.begin(), minimizer_list.begin() + safe_ct);
      index_list.erase(index_list.begin(), index_list.begin() + safe_ct);
      insertion_list.erase(insertion_list.begin(), insertion_list.begin() + safe_ct);
      mm_ct -= safe_ct;
    }  // end deterministic order while loop
  }  // end block for loop

  for (int i = 0; i < set_ct; i++)
    omp_destroy_lock(&locks[i]);
}

#endif
