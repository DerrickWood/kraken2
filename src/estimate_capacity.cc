/*
 * Copyright 2013-2020, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include "kraken2_headers.h"
#include "kv_store.h"
#include "mmscanner.h"
#include "seqreader.h"
#include "utilities.h"

using std::string;
using std::cout;
using std::queue;
using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::unordered_set;
using std::vector;
using namespace kraken2;

#define RANGE_SECTIONS 1024  // must be power of 2
#define RANGE_MASK (RANGE_SECTIONS - 1)
#define MAX_N RANGE_SECTIONS
#define DEFAULT_N 4
#define DEFAULT_BLOCK_SIZE (30 * 1024 * 1024)  // yes, 30 MB

struct Options {
  size_t k, l, n;
  bool input_is_protein;
  int threads;
  size_t block_size;
  uint64_t spaced_seed_mask;
  uint64_t toggle_mask;
};

void ParseCommandLine(int argc, char **argv, Options &opts);
void usage(int exit_code = EX_USAGE);
void ProcessSequence(string &seq, Options &opts,
    vector<unordered_set<uint64_t>> &sets);
void ProcessSequences(Options &opts);

int main(int argc, char **argv) {
  Options opts;
  opts.k = 0;
  opts.l = 0;
  opts.n = DEFAULT_N;
  opts.threads = 1;
  opts.input_is_protein = false;
  opts.spaced_seed_mask = DEFAULT_SPACED_SEED_MASK;
  opts.toggle_mask = DEFAULT_TOGGLE_MASK;
  opts.block_size = DEFAULT_BLOCK_SIZE;
  ParseCommandLine(argc, argv, opts);
  omp_set_num_threads(opts.threads);
  ProcessSequences(opts);
  return 0;
}

void ProcessSequences(Options &opts)
{
  vector<unordered_set<uint64_t>> sets(opts.n);

  #pragma omp parallel
  {
    bool have_work = true;
    BatchSequenceReader reader;
    Sequence sequence;

    while (have_work) {
      #pragma omp critical(batch_reading)
      have_work = reader.LoadBlock(std::cin, opts.block_size);
      if (have_work)
        while (reader.NextSequence(sequence))
          ProcessSequence(sequence.seq, opts, sets);
    }
  }

  size_t sum_set_sizes = 0;
  for (auto &s : sets) {
    sum_set_sizes += s.size();
  }
  sum_set_sizes++;  // ensure non-zero estimate
  cout << (size_t) (sum_set_sizes * RANGE_SECTIONS * 1.0 / opts.n) << endl;
}

void ParseCommandLine(int argc, char **argv, Options &opts) {
  int opt;
  long long sig;

  while ((opt = getopt(argc, argv, "?hk:l:n:S:T:B:p:X")) != -1) {
    switch (opt) {
      case 'h' : case '?' :
        usage(0);
        break;
      case 'p' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "must have at least 1 thread");
        opts.threads = sig;
        break;
      case 'B' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "block size must be positive");
        opts.block_size = sig;
        break;
      case 'n' :
        sig = atoll(optarg);
        if (sig < 1)
          errx(EX_USAGE, "n must be positive integer");
        if (sig > MAX_N)
          errx(EX_USAGE, "n must be no more than %d", MAX_N);
        opts.n = sig;
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
      case 'X' :
        opts.input_is_protein = true;
        break;
      case 'S' :
        opts.spaced_seed_mask = strtol(optarg, nullptr, 2);
        break;
      case 'T' :
        opts.toggle_mask = strtol(optarg, nullptr, 2);
        break;
    }
  }

  if (opts.spaced_seed_mask != DEFAULT_SPACED_SEED_MASK)
    ExpandSpacedSeedMask(opts.spaced_seed_mask,
      opts.input_is_protein ? BITS_PER_CHAR_PRO : BITS_PER_CHAR_DNA);
  if (opts.k == 0 || opts.l == 0) {
    cerr << "missing mandatory integer parameter" << endl;
    usage();
  }
  if (opts.k < opts.l) {
    cerr << "k cannot be less than l" << endl;
    usage();
  }
}

void usage(int exit_code) {
  cerr << "Usage: estimate_capacity <options>" << endl
       << endl
       << "Options (*mandatory):" << endl
       << "* -k INT        Set length of k-mers" << endl
       << "* -l INT        Set length of minimizers" << endl
       << "  -n INT        Set maximum qualifying hash code" << endl
       << "  -X            Input sequences are proteins" << endl
       << "  -S BITSTRING  Spaced seed mask" << endl
       << "  -T BITSTRING  Minimizer ordering toggle mask" << endl
       << "  -B INT        Read block size" << endl
       << "  -p INT        Number of threads" << endl;
  exit(exit_code);
}

void ProcessSequence(string &seq, Options &opts,
    vector<unordered_set<uint64_t>> &sets)
{
  MinimizerScanner scanner(opts.k, opts.l, opts.spaced_seed_mask, ! opts.input_is_protein, opts.toggle_mask);
  // Add terminator for protein sequences if not already there
  if (opts.input_is_protein && seq.back() != '*')
    seq.push_back('*');
  scanner.LoadSequence(seq);
  uint64_t *minimizer_ptr;
  while ((minimizer_ptr = scanner.NextMinimizer())) {
    if (scanner.is_ambiguous())
      continue;
    uint64_t hash_code = MurmurHash3(*minimizer_ptr);
    if ((hash_code & RANGE_MASK) < opts.n) {
      #pragma omp critical(set_insert)
      sets[hash_code & RANGE_MASK].insert(*minimizer_ptr);
    }
  }
}
