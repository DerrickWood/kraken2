/*
 * Copyright 2013-2020, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include "kraken2_headers.h"
#include "taxonomy.h"
#include "mmscanner.h"
#include "seqreader.h"
#include "compact_hash.h"
#include "kv_store.h"
#include "kraken2_data.h"
#include "utilities.h"

using std::string;
using std::map;
using std::cout;
using std::queue;
using std::cerr;
using std::endl;
using std::vector;
using std::istringstream;
using std::ifstream;
using std::ofstream;
using namespace kraken2;

#define DEFAULT_BLOCK_SIZE (10 * 1024 * 1024)  // 10 MB

struct Options {
  string ID_to_taxon_map_filename;
  string ncbi_taxonomy_directory;
  string hashtable_filename;
  string options_filename;
  string taxonomy_filename;
  size_t block_size;
  int num_threads;
  bool input_is_protein;
  ssize_t k, l;
  size_t capacity;
  size_t maximum_capacity;
  uint64_t spaced_seed_mask;
  uint64_t toggle_mask;
  uint64_t min_clear_hash_value;
};

void ParseCommandLine(int argc, char **argv, Options &opts);
void usage(int exit_code = EX_USAGE);
vector<string> ExtractNCBISequenceIDs(const string &header);
void ProcessSequence(string &seq, uint64_t taxid,
    CompactHashTable &hash, Taxonomy &tax, MinimizerScanner &scanner,
    uint64_t min_clear_hash_value);
void ProcessSequences(Options &opts, map<string, uint64_t> &ID_to_taxon_map,
    CompactHashTable &kraken_index, Taxonomy &taxonomy);
void ReadIDToTaxonMap(map<string, uint64_t> &id_map, string &filename);
void GenerateTaxonomy(Options &opts, map<string, uint64_t> &id_map);

struct TaxonSeqPair {
  uint64_t taxon;
  string seq;
};

int main(int argc, char **argv) {
  Options opts;
  opts.spaced_seed_mask = DEFAULT_SPACED_SEED_MASK;
  opts.toggle_mask = DEFAULT_TOGGLE_MASK;
  opts.input_is_protein = false;
  opts.num_threads = 1;
  opts.block_size = DEFAULT_BLOCK_SIZE;
  opts.min_clear_hash_value = 0;
  opts.maximum_capacity = 0;
  ParseCommandLine(argc, argv, opts);

  omp_set_num_threads( opts.num_threads );

  map<string, uint64_t> ID_to_taxon_map;

  ReadIDToTaxonMap(ID_to_taxon_map, opts.ID_to_taxon_map_filename);
  GenerateTaxonomy(opts, ID_to_taxon_map);

  std::cerr << "Taxonomy parsed and converted." << std::endl;

  Taxonomy taxonomy(opts.taxonomy_filename.c_str());
  taxonomy.GenerateExternalToInternalIDMap();
  size_t bits_needed_for_value = 1;
  while ((1 << bits_needed_for_value) < (ssize_t) taxonomy.node_count())
    bits_needed_for_value++;

  auto actual_capacity = opts.capacity;
  if (opts.maximum_capacity) {
    double frac = opts.maximum_capacity * 1.0 / opts.capacity;
    if (frac > 1)
      errx(EX_DATAERR, "maximum capacity larger than requested capacity");
    opts.min_clear_hash_value = (uint64_t) ((1 - frac) * UINT64_MAX);
    actual_capacity = opts.maximum_capacity;
  }

  CompactHashTable kraken_index(actual_capacity, 32 - bits_needed_for_value,
      bits_needed_for_value);
  std::cerr << "CHT created with " << bits_needed_for_value << " bits reserved for taxid." << std::endl;

  ProcessSequences(opts, ID_to_taxon_map, kraken_index, taxonomy);

  std::cerr << "Writing data to disk... " << std::flush;
  kraken_index.WriteTable(opts.hashtable_filename.c_str());

  IndexOptions index_opts;
  index_opts.k = opts.k;
  index_opts.l = opts.l;
  index_opts.spaced_seed_mask = opts.spaced_seed_mask;
  index_opts.toggle_mask = opts.toggle_mask;
  index_opts.dna_db = ! opts.input_is_protein;
  index_opts.minimum_acceptable_hash_value = opts.min_clear_hash_value;
  index_opts.revcom_version = CURRENT_REVCOM_VERSION;
  ofstream opts_fs(opts.options_filename);
  opts_fs.write((char *) &index_opts, sizeof(index_opts));
  if (! opts_fs.good())
    errx(EX_OSERR, "Unable to write options file %s",
         opts.options_filename.c_str());
  opts_fs.close();
  std::cerr << " complete." << std::endl;

  return 0;
}

void ProcessSequences(Options &opts, map<string, uint64_t> &ID_to_taxon_map,
    CompactHashTable &kraken_index, Taxonomy &taxonomy)
{
  size_t processed_seq_ct = 0;
  size_t processed_ch_ct = 0;

  #pragma omp parallel
  {
    Sequence sequence;
    MinimizerScanner scanner(opts.k, opts.l, opts.spaced_seed_mask,
                             ! opts.input_is_protein, opts.toggle_mask);

    BatchSequenceReader reader;

    while (true) {
      // Declaration of "ok" and break need to be done outside of critical
      // section to conform with OpenMP spec.
      bool ok;
      #pragma omp critical(reader)
      ok = reader.LoadBlock(std::cin, opts.block_size);
      if (! ok)
        break;
      while (reader.NextSequence(sequence)) {
        auto all_sequence_ids = ExtractNCBISequenceIDs(sequence.header);
        uint64_t taxid = 0;
        for (auto &seqid : all_sequence_ids) {
          if (ID_to_taxon_map.count(seqid) == 0) continue;
          auto ext_taxid = ID_to_taxon_map[seqid];
          taxid = taxonomy.LowestCommonAncestor(taxid, taxonomy.GetInternalID(ext_taxid));
        }
        if (taxid) {
          // Add terminator for protein sequences if not already there
          if (opts.input_is_protein && sequence.seq.back() != '*')
            sequence.seq.push_back('*');
          ProcessSequence(sequence.seq, taxid, kraken_index, taxonomy, scanner,
            opts.min_clear_hash_value);
          #pragma omp atomic
          processed_seq_ct++;
          #pragma omp atomic
          processed_ch_ct += sequence.seq.size();
        }
      }
      if (isatty(fileno(stderr))) {
        #pragma omp critical(status_update)
        std::cerr << "\rProcessed " << processed_seq_ct << " sequences (" << processed_ch_ct << " " << (opts.input_is_protein ? "aa" : "bp") << ")...";
      }
    }
  }
  if (isatty(fileno(stderr)))
    std::cerr << "\r";
  std::cerr << "Completed processing of " << processed_seq_ct << " sequences, " << processed_ch_ct << " " << (opts.input_is_protein ? "aa" : "bp") << std::endl;
}

// This function exists to deal with NCBI's use of \x01 characters to denote
// the start of a new FASTA header in the same line (for non-redundant DBs).
// We return all sequence IDs in a header line, not just the first.
vector<string> ExtractNCBISequenceIDs(const string &header) {
  vector<string> list;
  string current_str;

  bool in_id = true;
  // start loop at first char after '>'
  for (size_t i = 1; i < header.size(); ++i) {
    if (header[i] == 0x01) {
      // 0x01 starts new ID at next char
      if (! current_str.empty())
        list.push_back(current_str);
      current_str.clear();
      in_id = true;
    }
    else if (in_id && isspace(header[i])) {
      // spaces end ID
      if (! current_str.empty())
        list.push_back(current_str);
      current_str.clear();
      in_id = false;
    }
    else if (in_id) {
      // Build ID string char by char
      current_str.push_back(header[i]);
    }
  }
  if (! current_str.empty())
    list.push_back(current_str);
  return list;
}

void ParseCommandLine(int argc, char **argv, Options &opts) {
  int opt;
  long long sig;

  while ((opt = getopt(argc, argv, "?hB:c:H:m:n:o:t:k:l:s:S:T:p:M:X")) != -1) {
    switch (opt) {
      case 'h' : case '?' :
        usage(0);
        break;
      case 'B' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "can't have negative block size");
        opts.block_size = sig;
        break;
      case 'p' :
        opts.num_threads = atoll(optarg);
        if (opts.num_threads < 1)
          errx(EX_USAGE, "can't have negative number of threads");
        if (opts.num_threads > omp_get_max_threads())
          errx(EX_USAGE, "OMP only wants you to use %d threads", omp_get_max_threads());
        break;
      case 'H' :
        opts.hashtable_filename = optarg;
        break;
      case 'm' :
        opts.ID_to_taxon_map_filename = optarg;
        break;
      case 'n' :
        opts.ncbi_taxonomy_directory = optarg;
        break;
      case 'o' :
        opts.options_filename = optarg;
        break;
      case 't' :
        opts.taxonomy_filename = optarg;
        break;
      case 'S' :
        opts.spaced_seed_mask = strtol(optarg, nullptr, 2);
        break;
      case 'T' :
        opts.toggle_mask = strtol(optarg, nullptr, 2);
        break;
      case 'k' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "k must be positive integer");
        opts.k = sig;
        break;
      case 'l' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "l must be positive integer");
        if (sig > 31)
          errx(EX_USAGE, "l must be no more than 31");
        opts.l = sig;
        break;
      case 'c' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "capacity must be positive integer");
        opts.capacity = sig;
        break;
      case 'M' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "max capacity must be positive integer");
        opts.maximum_capacity = sig;
        break;
      case 'X' :
        opts.input_is_protein = true;
        break;
    }
  }

  if (opts.spaced_seed_mask != DEFAULT_SPACED_SEED_MASK)
    ExpandSpacedSeedMask(opts.spaced_seed_mask,
      opts.input_is_protein ? BITS_PER_CHAR_PRO : BITS_PER_CHAR_DNA);
  if (opts.hashtable_filename.empty() ||
      opts.ID_to_taxon_map_filename.empty() ||
      opts.ncbi_taxonomy_directory.empty() ||
      opts.options_filename.empty() ||
      opts.taxonomy_filename.empty())
  {
    cerr << "missing mandatory filename parameter" << endl;
    usage();
  }
  if (opts.k == 0 || opts.l == 0 || opts.capacity == 0) {
    cerr << "missing mandatory integer parameter" << endl;
    usage();
  }
  if (opts.k < opts.l) {
    cerr << "k cannot be less than l" << endl;
    usage();
  }
  if (opts.maximum_capacity > opts.capacity) {
    cerr << "maximum capacity option shouldn't specify larger capacity than normal" << endl;
    usage();
  }
}

void usage(int exit_code) {
  cerr << "Usage: build_db <options>\n"
       << "\n"
       << "Options (*mandatory):\n"
       << "* -H FILENAME   Kraken 2 hash table filename\n"
       << "* -m FILENAME   Sequence ID to taxon map filename\n"
       << "* -t FILENAME   Kraken 2 taxonomy filename\n"
       << "* -n DIR        NCBI taxonomy directory name\n"
       << "* -o FILENAME   Kraken 2 options filename\n"
       << "* -k INT        Set length of k-mers\n"
       << "* -l INT        Set length of minimizers\n"
       << "* -c INT        Set capacity of hash table\n"
       << "  -M INT        Set maximum capacity of hash table (MiniKraken)\n"
       << "  -S BITSTRING  Spaced seed mask\n"
       << "  -T BITSTRING  Minimizer toggle mask\n"
       << "  -X            Input seqs. are proteins\n"
       << "  -p INT        Number of threads\n"
       << "  -B INT        Read block size" << endl;
  exit(exit_code);
}

void ProcessSequence(string &seq, uint64_t taxid,
    CompactHashTable &hash, Taxonomy &tax, MinimizerScanner &scanner,
    uint64_t min_clear_hash_value)
{
  scanner.LoadSequence(seq);
  uint64_t *minimizer_ptr;
  while ((minimizer_ptr = scanner.NextMinimizer())) {
    if (min_clear_hash_value && MurmurHash3(*minimizer_ptr) < min_clear_hash_value)
      continue;
    hvalue_t existing_taxid = 0;
    hvalue_t new_taxid = taxid;
    while (! hash.CompareAndSet(*minimizer_ptr, new_taxid, &existing_taxid)) {
      new_taxid = tax.LowestCommonAncestor(new_taxid, existing_taxid);
    }
  }
}

void ReadIDToTaxonMap(map<string, uint64_t> &id_map, string &filename) {
  ifstream map_file(filename.c_str());
  if (! map_file.good())
    err(EX_NOINPUT, "unable to read from '%s'", filename.c_str());
  string line;

  while (getline(map_file, line)) {
    string seq_id;
    uint64_t taxid;
    istringstream iss(line);
    iss >> seq_id;
    iss >> taxid;
    id_map[seq_id] = taxid;
  }
}

void GenerateTaxonomy(Options &opts, map<string, uint64_t> &id_map) {
  NCBITaxonomy ncbi_taxonomy(
    opts.ncbi_taxonomy_directory + "/nodes.dmp",
    opts.ncbi_taxonomy_directory + "/names.dmp"
  );
  for (auto &kv_pair : id_map) {
    if (kv_pair.second != 0) {
      ncbi_taxonomy.MarkNode(kv_pair.second);
    }
  }
  ncbi_taxonomy.ConvertToKrakenTaxonomy(opts.taxonomy_filename.c_str());
}
