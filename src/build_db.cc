/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include "build_db.h"
#include "kraken2_headers.h"
#include "taxonomy.h"
#include "mmscanner.h"
#include "compact_hash.h"
#include "kraken2_data.h"
#include "utilities.h"

using std::string;
using std::map;
using std::cerr;
using std::endl;
using std::vector;
using std::istringstream;
using std::ifstream;
using std::ofstream;
using namespace kraken2;

// These will likely be overridden by wrapper script,
// just keeping sane defaults in case they aren't
#define DEFAULT_BLOCK_SIZE (10 * 1024 * 1024)  // 10 MB
#define DEFAULT_SUBBLOCK_SIZE (1024)

void ParseCommandLine(int argc, char **argv, Options &opts);
void usage(int exit_code = EX_USAGE);
void ReadIDToTaxonMap(map<string, taxid_t> &id_map, string &filename);
void GenerateTaxonomy(Options &opts, map<string, taxid_t> &id_map);

struct TaxonSeqPair {
  taxid_t taxon;
  string seq;
};

int main(int argc, char **argv) {
  Options opts;
  opts.spaced_seed_mask = DEFAULT_SPACED_SEED_MASK;
  opts.toggle_mask = DEFAULT_TOGGLE_MASK;
  opts.input_is_protein = false;
  opts.num_threads = 1;
  opts.block_size = DEFAULT_BLOCK_SIZE;
  opts.subblock_size = DEFAULT_SUBBLOCK_SIZE;
  opts.cht_cell_size = 32;
  opts.requested_bits_for_taxid = 0;
  opts.min_clear_hash_value = 0;
  opts.maximum_capacity = 0;
  opts.deterministic_build = true;
  ParseCommandLine(argc, argv, opts);

  omp_set_num_threads( opts.num_threads );

  map<string, taxid_t> ID_to_taxon_map;

  ReadIDToTaxonMap(ID_to_taxon_map, opts.ID_to_taxon_map_filename);
  GenerateTaxonomy(opts, ID_to_taxon_map);

  std::cerr << "Taxonomy parsed and converted." << std::endl;

  Taxonomy taxonomy(opts.taxonomy_filename.c_str(), false);
  taxonomy.GenerateExternalToInternalIDMap();
  size_t bits_needed_for_value = 1;
  while ((1 << bits_needed_for_value) < (ssize_t) taxonomy.node_count())
    bits_needed_for_value++;
  if (opts.requested_bits_for_taxid > 0 &&
      bits_needed_for_value > opts.requested_bits_for_taxid)
    errx(EX_DATAERR, "more bits required for storing taxid");

  size_t bits_for_taxid = bits_needed_for_value;
  if (bits_for_taxid < opts.requested_bits_for_taxid)
    bits_for_taxid = opts.requested_bits_for_taxid;

  auto actual_capacity = opts.capacity;
  if (opts.maximum_capacity) {
    double frac = opts.maximum_capacity * 1.0 / opts.capacity;
    if (frac > 1)
      errx(EX_DATAERR, "maximum capacity larger than requested capacity");
    opts.min_clear_hash_value = (uint64_t) ((1 - frac) * UINT64_MAX);
    actual_capacity = opts.maximum_capacity;
  }

  if (opts.cht_cell_size == 32) {
    build<CompactHashCell40>(taxonomy, ID_to_taxon_map, opts, actual_capacity, bits_for_taxid);
  } else if (opts.cht_cell_size == 40) {
    build<CompactHashCell40>(taxonomy, ID_to_taxon_map, opts, actual_capacity, bits_for_taxid);
  } else {
    errx(EX_DATAERR, "Unsupported CHT cell size");
  }

  IndexOptions index_opts;
  index_opts.k = opts.k;
  index_opts.l = opts.l;
  index_opts.spaced_seed_mask = opts.spaced_seed_mask;
  index_opts.toggle_mask = opts.toggle_mask;
  index_opts.dna_db = ! opts.input_is_protein;
  index_opts.minimum_acceptable_hash_value = opts.min_clear_hash_value;
  index_opts.revcom_version = CURRENT_REVCOM_VERSION;
  index_opts.db_version = 0;
  index_opts.db_type = 0;
  ofstream opts_fs(opts.options_filename);
  opts_fs.write((char *) &index_opts, sizeof(index_opts));
  if (! opts_fs.good())
    errx(EX_OSERR, "Unable to write options file %s",
         opts.options_filename.c_str());
  opts_fs.close();
  std::cerr << " complete." << std::endl;

  return 0;
}

// This function exists to deal with NCBI's use of \x01 characters to denote
// the start of a new FASTA header in the same line (for non-redundant DBs).
// We return all sequence IDs in a header line, not just the first.
vector<string> ExtractNCBISequenceIDs(const string &header) {
  vector<string> list;
  string current_str;
  bool in_id = true;

  for (size_t i = 0; i < header.size(); ++i) {
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

void ReadIDToTaxonMap(map<string, taxid_t> &id_map, string &filename) {
  ifstream map_file(filename.c_str());
  if (! map_file.good())
    err(EX_NOINPUT, "unable to read from '%s'", filename.c_str());
  string line;

  while (getline(map_file, line)) {
    string seq_id;
    taxid_t taxid;
    istringstream iss(line);
    iss >> seq_id;
    iss >> taxid;
    if (taxid)
      id_map[seq_id] = taxid;
  }
}

void GenerateTaxonomy(Options &opts, map<string, taxid_t> &id_map) {
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

void ParseCommandLine(int argc, char **argv, Options &opts) {
  int opt;
  long long sig;

  while ((opt = getopt(argc, argv, "?hB:b:c:FH:m:n:o:t:k:l:M:p:r:s:S:T:X")) != -1) {
    switch (opt) {
      case 'h' : case '?' :
        usage(0);
        break;
      case 'B' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "must have positive block size");
        opts.block_size = sig;
        break;
      case 'b' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "must have positive subblock size");
        opts.subblock_size = sig;
        break;
      case 'r' :
        sig = atoll(optarg);
        if (sig < 0)
          errx(EX_USAGE, "can't have negative bit storage");
        if (sig > 31)
          errx(EX_USAGE, "can't have more than 31 bits of storage for taxid");
        opts.requested_bits_for_taxid = sig;
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
      case 'F' :
        opts.deterministic_build = false;
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
  if (opts.block_size < opts.subblock_size) {
    cerr << "block size cannot be less than subblock size\n";
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
       << "  -F            Use fast, nondeterministic building method\n"
       << "  -B INT        Read block size\n"
       << "  -b INT        Read subblock size\n"
       << "  -r INT        Bit storage requested for taxid" << endl;
  exit(exit_code);
}
