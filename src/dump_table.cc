/*
 * Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include "kraken2_headers.h"
#include "compact_hash.h"
#include "taxonomy.h"
#include "mmscanner.h"
#include "kraken2_data.h"
#include "reports.h"

using std::string;
using std::unordered_map;
using std::cout;
using std::cerr;
using std::endl;
using namespace kraken2;

struct Options {
  string hashtable_filename;
  string taxonomy_filename;
  string options_filename;
  string output_filename;
  bool use_mpa_style;
  bool report_zeros;
  bool skip_counts;
  unsigned int num_threads;
};

void ParseCommandLine(int argc, char **argv, Options &opts);
void usage(int exit_code = EX_USAGE);

std::string mask2str(uint64_t mask, int digits) {
  std::string str;
  for (int i = digits - 1; i >= 0; i--)
    str.push_back((mask >> i) & 1 ? '1' : '0');
  return str;
}

int main(int argc, char **argv) {
  Options opts;
  opts.report_zeros = false;
  opts.output_filename = "/dev/fd/1";
  opts.use_mpa_style = false;
  opts.skip_counts = false;
  opts.num_threads = 1;
  ParseCommandLine(argc, argv, opts);

  omp_set_num_threads(opts.num_threads);

  CompactHashTable kraken_index(opts.hashtable_filename);
  Taxonomy taxonomy(opts.taxonomy_filename);
  IndexOptions idx_opts = {0};
  std::ifstream idx_opt_fs(opts.options_filename);
  idx_opt_fs.read((char *) &idx_opts, sizeof(idx_opts));

  std::cout << "# Database options: "
    << (idx_opts.dna_db ? "nucleotide" : "protein") << " db, "
    << "k = " << idx_opts.k << ", "
    << "l = " << idx_opts.l << "\n"
    << "# Spaced mask = "
      << mask2str(idx_opts.spaced_seed_mask,
          idx_opts.l * (idx_opts.dna_db ? BITS_PER_CHAR_DNA : BITS_PER_CHAR_PRO))
      << "\n"
    << "# Toggle mask = " << mask2str(idx_opts.toggle_mask, 64) << "\n"
    << "# Total taxonomy nodes: " << taxonomy.node_count() << "\n"
    << "# Table size: " << kraken_index.size() << "\n"
    << "# Table capacity: " << kraken_index.capacity() << "\n"
    << "# Min clear hash value = " << idx_opts.minimum_acceptable_hash_value << "\n";
  if (idx_opts.revcom_version != CURRENT_REVCOM_VERSION)
    std::cout << "# Built with outdated revcom version\n";
  std::cout << std::flush;

  if (! opts.skip_counts) {
    taxon_counts_t taxid_counts = kraken_index.GetValueCounts();
    uint64_t total_seqs = 0;
    taxon_counters_t taxid_counters;
    for (auto &kv_pair : taxid_counts) {
      total_seqs += kv_pair.second;
      READCOUNTER rc(kv_pair.second, 0);
      taxid_counters[kv_pair.first] = rc;
    }
    if (opts.use_mpa_style)
      ReportMpaStyle(opts.output_filename, opts.report_zeros, taxonomy, taxid_counters);
    else
      ReportKrakenStyle(opts.output_filename, opts.report_zeros, false, taxonomy,
          taxid_counters, total_seqs, 0);
  }

  return 0;
}

void ParseCommandLine(int argc, char **argv, Options &opts) {
  int opt;

  while ((opt = getopt(argc, argv, "?hH:t:o:O:p:zms")) != -1) {
    switch (opt) {
      case 'h' : case '?' :
        usage(0);
        break;
      case 'H' :
        opts.hashtable_filename = optarg;
        break;
      case 't' :
        opts.taxonomy_filename = optarg;
        break;
      case 'o' :
        opts.options_filename = optarg;
        break;
      case 'z' :
        opts.report_zeros = true;
        break;
      case 'O' :
        opts.output_filename = optarg;
        break;
      case 'm' :
        opts.use_mpa_style = true;
        break;
      case 's' :
        opts.skip_counts = true;
        break;
      case 'p' :
        opts.num_threads = atoi(optarg);
        break;
    }
  }

  if (opts.hashtable_filename.empty() || opts.taxonomy_filename.empty() ||
      opts.options_filename.empty())
  {
    cerr << "missing mandatory filename parameter" << endl;
    usage();
  }
}

void usage(int exit_code) {
  cerr << "Usage: dump_table <options>\n"
       << "\n"
       << "Options (*mandatory):\n"
       << "* -H FILENAME   Kraken 2 hash table filename\n"
       << "* -t FILENAME   Kraken 2 taxonomy filename\n"
       << "* -o FILENAME   Kraken 2 database options filename\n"
       << "  -O FILENAME   Output filename (def: /dev/fd/1)\n"
       << "  -m            Use MPA style output instead of Kraken 2 output\n"
       << "  -s            Skip reporting minimizer counts, just show DB stats\n"
       << "  -z            Report taxa with zero counts\n";
  exit(exit_code);
}
