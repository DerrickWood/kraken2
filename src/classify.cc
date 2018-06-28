/*
 * Copyright 2013-2018, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include "kraken2_headers.h"
#include "kv_store.h"
#include "taxonomy.h"
#include "seqreader.h"
#include "mmscanner.h"
#include "compact_hash.h"
#include "kraken2_data.h"
#include "aa_translate.h"
#include "reports.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::map;
using std::ostringstream;
using std::ofstream;
using std::set;
using std::string;
using std::vector;
using namespace kraken2;

static const size_t NUM_FRAGMENTS_PER_THREAD = 10000;
static const taxid_t MATE_PAIR_BORDER_TAXON = TAXID_MAX;
static const taxid_t READING_FRAME_BORDER_TAXON = TAXID_MAX - 1;
static const taxid_t AMBIGUOUS_SPAN_TAXON = TAXID_MAX - 2;

struct Options {
  string index_filename;
  string taxonomy_filename;
  string options_filename;
  string report_filename;
  string classified_output_prefix;
  string unclassified_output_prefix;
  string kraken_output_filename;
  bool mpa_style_report;
  bool quick_mode;
  bool report_zero_counts;
  bool use_translated_search;
  bool print_scientific_name;
  double confidence_threshold;
  int num_threads;
  bool paired_end_processing;
  bool single_file_pairs;
  int minimum_quality_score;
  bool use_memory_mapping;
};

struct ClassificationStats {
  uint64_t total_sequences;
  uint64_t total_bases;
  uint64_t total_classified;
};

struct OutputStreamData {
  bool initialized;
  bool printing_sequences;
  std::ostream *classified_output1;
  std::ostream *classified_output2;
  std::ostream *unclassified_output1;
  std::ostream *unclassified_output2;
  std::ostream *kraken_output;
};

struct OutputData {
  uint64_t block_id;
  string kraken_str;
  string classified_out1_str;
  string classified_out2_str;
  string unclassified_out1_str;
  string unclassified_out2_str;
};

void ParseCommandLine(int argc, char **argv, Options &opts);
void usage(int exit_code=EX_USAGE);
void ProcessFiles(const char *filename1, const char *filename2,
    KeyValueStore *hash, Taxonomy &tax,
    IndexOptions &idx_opts, Options &opts, ClassificationStats &stats,
    OutputStreamData &outputs, taxon_counts_t &call_counts);
taxid_t ClassifySequence(Sequence &dna, Sequence &dna2, ostringstream &koss,
    KeyValueStore *hash, Taxonomy &tax, IndexOptions &idx_opts,
    Options &opts, ClassificationStats &stats, MinimizerScanner &scanner,
    vector<taxid_t> &taxa, taxon_counts_t &hit_counts,
    vector<string> &tx_frames);
void AddHitlistString(ostringstream &oss, vector<taxid_t> &taxa,
    Taxonomy &taxonomy);
taxid_t ResolveTree(taxon_counts_t &hit_counts,
    Taxonomy &tax, size_t total_minimizers, Options &opts);
void ReportStats(struct timeval time1, struct timeval time2,
    ClassificationStats &stats);
void InitializeOutputs(Options &opts, OutputStreamData &outputs, SequenceFormat format);
void MaskLowQualityBases(Sequence &dna, int minimum_quality_score);

int main(int argc, char **argv) {
  Options opts;
  opts.quick_mode = false;
  opts.confidence_threshold = 0;
  opts.paired_end_processing = false;
  opts.single_file_pairs = false;
  opts.num_threads = 1;
  opts.mpa_style_report = false;
  opts.report_zero_counts = false;
  opts.use_translated_search = false;
  opts.print_scientific_name = false;
  opts.minimum_quality_score = 0;
  opts.use_memory_mapping = false;

  ParseCommandLine(argc, argv, opts);

  omp_set_num_threads(opts.num_threads);

  cerr << "Loading database information...";

  IndexOptions idx_opts;
  ifstream idx_opt_fs(opts.options_filename);
  idx_opt_fs.read((char *) &idx_opts, sizeof(idx_opts));
  opts.use_translated_search = ! idx_opts.dna_db;

  Taxonomy taxonomy(opts.taxonomy_filename, opts.use_memory_mapping);
  KeyValueStore *hash_ptr = new CompactHashTable(opts.index_filename, opts.use_memory_mapping);

  cerr << " done." << endl;

  ClassificationStats stats = {0, 0, 0};
  taxon_counts_t call_counts;

  OutputStreamData outputs = { false, false, nullptr, nullptr, nullptr, nullptr, &std::cout };

  struct timeval tv1, tv2;
  gettimeofday(&tv1, nullptr);
  if (optind == argc) {
    if (opts.paired_end_processing && ! opts.single_file_pairs)
      errx(EX_USAGE, "paired end processing used with no files specified");
    ProcessFiles(nullptr, nullptr, hash_ptr, taxonomy, idx_opts, opts, stats, outputs, call_counts);
  }
  else {
    for (int i = optind; i < argc; i++) {
      if (opts.paired_end_processing && ! opts.single_file_pairs) {
        if (i + 1 == argc) {
          errx(EX_USAGE, "paired end processing used with unpaired file");
        }
        ProcessFiles(argv[i], argv[i+1], hash_ptr, taxonomy, idx_opts, opts, stats, outputs, call_counts);
        i += 1;
      }
      else {
        ProcessFiles(argv[i], nullptr, hash_ptr, taxonomy, idx_opts, opts, stats, outputs, call_counts);
      }
    }
  }
  gettimeofday(&tv2, nullptr);

  delete hash_ptr;

  ReportStats(tv1, tv2, stats);

  if (! opts.report_filename.empty()) {
    if (opts.mpa_style_report)
      ReportMpaStyle(opts.report_filename, opts.report_zero_counts, taxonomy,
          call_counts);
    else {
      auto total_unclassified = stats.total_sequences - stats.total_classified;
      ReportKrakenStyle(opts.report_filename, opts.report_zero_counts, taxonomy,
          call_counts, stats.total_sequences, total_unclassified);
    }
  }

  return 0;
}

void ReportStats(struct timeval time1, struct timeval time2,
  ClassificationStats &stats)
{
  time2.tv_usec -= time1.tv_usec;
  time2.tv_sec -= time1.tv_sec;
  if (time2.tv_usec < 0) {
    time2.tv_sec--;
    time2.tv_usec += 1000000;
  }
  double seconds = time2.tv_usec;
  seconds /= 1e6;
  seconds += time2.tv_sec;

  uint64_t total_unclassified = stats.total_sequences - stats.total_classified;

  cerr << "\r";
  fprintf(stderr,
          "%llu sequences (%.2f Mbp) processed in %.3fs (%.1f Kseq/m, %.2f Mbp/m).\n",
          (unsigned long long) stats.total_sequences,
          stats.total_bases / 1.0e6,
          seconds,
          stats.total_sequences / 1.0e3 / (seconds / 60),
          stats.total_bases / 1.0e6 / (seconds / 60) );
  fprintf(stderr, "  %llu sequences classified (%.2f%%)\n",
          (unsigned long long) stats.total_classified,
          stats.total_classified * 100.0 / stats.total_sequences);
  fprintf(stderr, "  %llu sequences unclassified (%.2f%%)\n",
          (unsigned long long) total_unclassified,
          total_unclassified * 100.0 / stats.total_sequences);
}

void ProcessFiles(const char *filename1, const char *filename2,
    KeyValueStore *hash, Taxonomy &tax,
    IndexOptions &idx_opts, Options &opts, ClassificationStats &stats,
    OutputStreamData &outputs, taxon_counts_t &call_counts)
{
  std::istream *fptr1 = nullptr, *fptr2 = nullptr;

  if (filename1 == nullptr)
    fptr1 = &std::cin;
  else {
    fptr1 = new std::ifstream(filename1);
  }
  if (opts.paired_end_processing && ! opts.single_file_pairs) {
    fptr2 = new std::ifstream(filename2);
  }

  // The priority queue for output is designed to ensure fragment data
  // is output in the same order it was input
  auto comparator = [](const OutputData &a, const OutputData &b) {
    return a.block_id > b.block_id;
  };
  std::priority_queue<OutputData, vector<OutputData>, decltype(comparator)>
    output_queue(comparator);
  uint64_t next_input_block_id = 0;
  uint64_t next_output_block_id = 0;
  omp_lock_t output_lock;
  omp_init_lock(&output_lock);

  #pragma omp parallel
  {
    MinimizerScanner scanner(idx_opts.k, idx_opts.l, idx_opts.toggle_mask,
                             idx_opts.dna_db, idx_opts.spaced_seed_mask);
    vector<taxid_t> taxa;
    taxon_counts_t hit_counts;
    ostringstream kraken_oss, c1_oss, c2_oss, u1_oss, u2_oss;
    ClassificationStats thread_stats = {0, 0, 0};
    vector<taxid_t> calls;
    vector<string> translated_frames(6);
    BatchSequenceReader reader1, reader2;
    Sequence seq1, seq2;
    uint64_t block_id;
    OutputData out_data;

    while (true) {
      thread_stats.total_sequences = 0;
      thread_stats.total_bases = 0;
      thread_stats.total_classified = 0;

      auto ok_read = false;

      #pragma omp critical(seqread)
      {  // Input processing block
        if (! opts.paired_end_processing) {
          // Unpaired data?  Just read in a sized block
          ok_read = reader1.LoadBlock(*fptr1, (size_t)(3 * 1024 * 1024));
        }
        else if (! opts.single_file_pairs) {
          // Paired data in 2 files?  Read a line-counted batch from each file.
          ok_read = reader1.LoadBatch(*fptr1, NUM_FRAGMENTS_PER_THREAD);
          if (ok_read && opts.paired_end_processing)
            ok_read = reader2.LoadBatch(*fptr2, NUM_FRAGMENTS_PER_THREAD);
        }
        else {
          auto frags = NUM_FRAGMENTS_PER_THREAD * 2;
          // Ensure frag count is even - just in case above line is changed
          if (frags % 2 == 1)
            frags++;
          ok_read = reader1.LoadBatch(*fptr1, frags);
        }
        block_id = next_input_block_id++;
      }

      if (! ok_read)
        break;

      // Reset all dynamically-growing things
      calls.clear();
      kraken_oss.str("");
      c1_oss.str("");
      c2_oss.str("");
      u1_oss.str("");
      u2_oss.str("");

      while (true) {
        auto valid_fragment = reader1.NextSequence(seq1);
        if (opts.paired_end_processing && valid_fragment) {
          if (opts.single_file_pairs)
            valid_fragment = reader1.NextSequence(seq2);
          else
            valid_fragment = reader2.NextSequence(seq2);
        }
        if (! valid_fragment)
          break;
        thread_stats.total_sequences++;
        if (opts.minimum_quality_score > 0) {
          MaskLowQualityBases(seq1, opts.minimum_quality_score);
          if (opts.paired_end_processing)
            MaskLowQualityBases(seq2, opts.minimum_quality_score);
        }
        auto call = ClassifySequence(seq1, seq2,
            kraken_oss, hash, tax, idx_opts, opts, thread_stats, scanner,
            taxa, hit_counts, translated_frames);
        if (call) {
          char buffer[1024] = "";
          sprintf(buffer, " kraken:taxid|%lu", call);
          seq1.header += buffer;
          seq2.header += buffer;
          c1_oss << seq1.to_string();
          if (opts.paired_end_processing)
            c2_oss << seq2.to_string();
        }
        else {
          u1_oss << seq1.to_string();
          if (opts.paired_end_processing)
            u2_oss << seq2.to_string();
        }
        calls.push_back(call);
        thread_stats.total_bases += seq1.seq.size();
        if (opts.paired_end_processing)
          thread_stats.total_bases += seq2.seq.size();
      }

      #pragma omp atomic
      stats.total_sequences += thread_stats.total_sequences;
      #pragma omp atomic
      stats.total_bases += thread_stats.total_bases;
      #pragma omp atomic
      stats.total_classified += thread_stats.total_classified;

      #pragma omp critical(output_stats)
      {
        cerr << "\rProcessed " << stats.total_sequences
             << " sequences (" << stats.total_bases << " bp) ...";
      }

      #pragma omp critical(update_calls)
      {
        for (auto &call : calls) {
          call_counts[call]++;
        }
      }

      if (! outputs.initialized) {
        InitializeOutputs(opts, outputs, reader1.file_format());
      }

      out_data.block_id = block_id;
      out_data.kraken_str.assign(kraken_oss.str());
      out_data.classified_out1_str.assign(c1_oss.str());
      out_data.classified_out2_str.assign(c2_oss.str());
      out_data.unclassified_out1_str.assign(u1_oss.str());
      out_data.unclassified_out2_str.assign(u2_oss.str());

      #pragma omp critical(output_queue)
      {
        output_queue.push(out_data);
      }

      bool output_loop = true;
      while (output_loop) {
        #pragma omp critical(output_queue)
        {
          output_loop = ! output_queue.empty();
          if (output_loop) {
            out_data = output_queue.top();
            if (out_data.block_id == next_output_block_id) {
              output_queue.pop();
              // Acquiring output lock obligates thread to print out
              // next output data block, contained in out_data
              omp_set_lock(&output_lock);
              next_output_block_id++;
            }
            else
              output_loop = false;
          }
        }
        if (! output_loop)
          break;
        // Past this point in loop, we know lock is set

        if (outputs.kraken_output != nullptr)
          (*outputs.kraken_output) << out_data.kraken_str;
        if (outputs.classified_output1 != nullptr)
          (*outputs.classified_output1) << out_data.classified_out1_str;
        if (outputs.classified_output2 != nullptr)
          (*outputs.classified_output2) << out_data.classified_out2_str;
        if (outputs.unclassified_output1 != nullptr)
          (*outputs.unclassified_output1) << out_data.unclassified_out1_str;
        if (outputs.unclassified_output2 != nullptr)
          (*outputs.unclassified_output2) << out_data.unclassified_out2_str;
        omp_unset_lock(&output_lock);
      }  // end while output loop
    }  // end while
  }  // end parallel block
  omp_destroy_lock(&output_lock);
  if (fptr1 != nullptr)
    delete fptr1;
  if (fptr2 != nullptr)
    delete fptr2;
  if (outputs.kraken_output != nullptr)
    (*outputs.kraken_output) << std::flush;
  if (outputs.classified_output1 != nullptr)
    (*outputs.classified_output1) << std::flush;
  if (outputs.classified_output2 != nullptr)
    (*outputs.classified_output2) << std::flush;
  if (outputs.unclassified_output1 != nullptr)
    (*outputs.unclassified_output1) << std::flush;
  if (outputs.unclassified_output2 != nullptr)
    (*outputs.unclassified_output2) << std::flush;
}

taxid_t ResolveTree(taxon_counts_t &hit_counts,
    Taxonomy &taxonomy, size_t total_minimizers, Options &opts)
{
  taxid_t max_taxon = 0;
  uint32_t max_score = 0;
  uint32_t required_score = ceil(opts.confidence_threshold * total_minimizers);

  // Sum each taxon's LTR path
  for (auto &kv_pair : hit_counts) {
    taxid_t taxon = kv_pair.first;
    uint32_t score = 0;

    for (auto &kv_pair2 : hit_counts) {
      taxid_t taxon2 = kv_pair2.first;

      if (taxonomy.IsAAncestorOfB(taxon2, taxon)) {
        score += kv_pair2.second;
      }
    }

    if (score > max_score) {
      max_score = score;
      max_taxon = taxon;
    }
    else if (score == max_score) {
      max_taxon = taxonomy.LowestCommonAncestor(max_taxon, taxon);
    }
  }

  // We probably have a call w/o required support (unless LCA resolved tie)
  while (max_taxon && max_score < required_score) {
    uint32_t score = 0;
    for (auto &kv_pair : hit_counts) {
      taxid_t taxon = kv_pair.first;
      // Add to score if taxon in max_taxon's clade
      if (taxonomy.IsAAncestorOfB(max_taxon, taxon)) {
        score += kv_pair.second;
      }
    }
    if (score >= required_score)
      // Kill loop and return, we've got enough support here
      return max_taxon;
    else
      // Run up tree until confidence threshold is met
      // Run off tree if required score isn't met
      max_taxon = taxonomy.nodes()[max_taxon].parent_id;
  }

  return max_taxon;
}

std::string TrimPairInfo(std::string &id) {
  size_t sz = id.size();
  if (sz <= 2)
    return id;
  if ( (id[sz - 2] == '_' || ispunct(id[sz - 2])) && isdigit(id[sz - 1]) )
    return id.substr(0, sz - 2);
  return id;
}

taxid_t ClassifySequence(Sequence &dna, Sequence &dna2, ostringstream &koss,
    KeyValueStore *hash, Taxonomy &taxonomy, IndexOptions &idx_opts,
    Options &opts, ClassificationStats &stats, MinimizerScanner &scanner,
    vector<taxid_t> &taxa, taxon_counts_t &hit_counts,
    vector<string> &tx_frames)
{
  uint64_t *minimizer_ptr;
  taxid_t call = 0;
  taxa.clear();
  hit_counts.clear();
  auto frame_ct = opts.use_translated_search ? 6 : 1;

  for (int mate_num = 0; mate_num < 2; mate_num++) {
    if (mate_num == 1 && ! opts.paired_end_processing)
      break;

    if (opts.use_translated_search) {
      TranslateToAllFrames(mate_num == 0 ? dna.seq : dna2.seq, tx_frames);
    }
    // index of frame is 0 - 5 w/ tx search (or 0 if no tx search)
    for (int frame_idx = 0; frame_idx < frame_ct; frame_idx++) {
      if (opts.use_translated_search) {
        scanner.LoadSequence(tx_frames[frame_idx]);
      }
      else {
        scanner.LoadSequence(mate_num == 0 ? dna.seq : dna2.seq);
      }
      uint64_t last_minimizer = UINT64_MAX;
      taxid_t last_taxon = TAXID_MAX;
      bool ambig_flag = false;
      while ((minimizer_ptr = scanner.NextMinimizer(&ambig_flag)) != nullptr) {
        taxid_t taxon;
        if (ambig_flag) {
          taxon = AMBIGUOUS_SPAN_TAXON;
        }
        else {
          if (*minimizer_ptr != last_minimizer) {
            taxon = hash->Get(*minimizer_ptr);
            last_taxon = taxon;
            last_minimizer = *minimizer_ptr;
          }
          else {
            taxon = last_taxon;
          }
          if (taxon) {
            if (opts.quick_mode) {
              call = taxon;
              goto finished_searching;  // need to break 3 loops here
            }
            hit_counts[taxon]++;
          }
        }
        taxa.push_back(taxon);
      }
      if (opts.use_translated_search && frame_idx != 5)
        taxa.push_back(READING_FRAME_BORDER_TAXON);
    }
    if (opts.paired_end_processing && mate_num == 0)
      taxa.push_back(MATE_PAIR_BORDER_TAXON);
  }

  finished_searching:

  auto total_kmers = taxa.size();
  if (opts.paired_end_processing)
    total_kmers--;  // account for the mate pair marker
  if (opts.use_translated_search)  // account for reading frame markers
    total_kmers -= opts.paired_end_processing ? 4 : 2;
  if (! opts.quick_mode)
    call = ResolveTree(hit_counts, taxonomy, total_kmers, opts);

  if (call)
    stats.total_classified++;

  if (call)
    koss << "C\t";
  else
    koss << "U\t";
  if (! opts.paired_end_processing)
    koss << dna.id << "\t";
  else
    koss << TrimPairInfo(dna.id) << "\t";

  auto ext_call = taxonomy.nodes()[call].external_id;
  if (opts.print_scientific_name) {
    const char *name = nullptr;
    if (call) {
      name = taxonomy.name_data() + taxonomy.nodes()[call].name_offset;
    }
    koss << (name ? name : "unclassified") << " (taxid " << ext_call << ")";
  }
  else {
    koss << ext_call;
  }

  koss << "\t";
  if (! opts.paired_end_processing)
    koss << dna.seq.size() << "\t";
  else
    koss << dna.seq.size() << "|" << dna2.seq.size() << "\t";

  if (opts.quick_mode) {
    koss << call << ":Q";
  }
  else {
    if (taxa.empty())
      koss << "0:0";
    else
      AddHitlistString(koss, taxa, taxonomy);
  }

  koss << endl;

  return call;
}

void AddHitlistString(ostringstream &oss, vector<taxid_t> &taxa,
    Taxonomy &taxonomy)
{
  auto last_code = taxa[0];
  auto code_count = 1;

  for (size_t i = 1; i < taxa.size(); i++) {
    auto code = taxa[i];

    if (code == last_code) {
      code_count += 1;
    }
    else {
      if (last_code != MATE_PAIR_BORDER_TAXON && last_code != READING_FRAME_BORDER_TAXON) {
        if (last_code == AMBIGUOUS_SPAN_TAXON) {
          oss << "A:" << code_count << " ";
        }
        else {
          auto ext_code = taxonomy.nodes()[last_code].external_id;
          oss << ext_code << ":" << code_count << " ";
        }
      }
      else {  // mate pair/reading frame marker
        oss << (last_code == MATE_PAIR_BORDER_TAXON ? "|:| " : "-:- ");
      }
      code_count = 1;
      last_code = code;
    }
  }
  if (last_code != MATE_PAIR_BORDER_TAXON && last_code != READING_FRAME_BORDER_TAXON) {
    if (last_code == AMBIGUOUS_SPAN_TAXON) {
      oss << "A:" << code_count << " ";
    }
    else {
      auto ext_code = taxonomy.nodes()[last_code].external_id;
      oss << ext_code << ":" << code_count;
    }
  }
  else {  // mate pair/reading frame marker
    oss << (last_code == MATE_PAIR_BORDER_TAXON ? "|:|" : "-:-");
  }
}

void InitializeOutputs(Options &opts, OutputStreamData &outputs, SequenceFormat format) {
  #pragma omp critical(output_init)
  {
    if (! outputs.initialized) {
      string extension = format == FORMAT_FASTA ? ".fa" : ".fq";
      if (! opts.classified_output_prefix.empty()) {
        string prefix = opts.classified_output_prefix;
        if (opts.paired_end_processing) {
          outputs.classified_output1 = new ofstream(prefix + "_1" + extension);
          outputs.classified_output2 = new ofstream(prefix + "_2" + extension);
        }
        else
          outputs.classified_output1 = new ofstream(prefix + extension);
        outputs.printing_sequences = true;
      }
      if (! opts.unclassified_output_prefix.empty()) {
        string prefix = opts.unclassified_output_prefix;
        if (opts.paired_end_processing) {
          outputs.unclassified_output1 = new ofstream(prefix + "_1" + extension);
          outputs.unclassified_output2 = new ofstream(prefix + "_2" + extension);
        }
        else
          outputs.unclassified_output1 = new ofstream(prefix + extension);
        outputs.printing_sequences = true;
      }
      if (! opts.kraken_output_filename.empty()) {
        if (opts.kraken_output_filename == "-")  // Special filename to silence Kraken output
          outputs.kraken_output = nullptr;
        else
          outputs.kraken_output = new ofstream(opts.kraken_output_filename);
      }
      outputs.initialized = true;
    }
  }
}

void MaskLowQualityBases(Sequence &dna, int minimum_quality_score) {
  if (dna.format != FORMAT_FASTQ)
    return;
  if (dna.seq.size() != dna.quals.size())
    errx(EX_DATAERR, "%s: Sequence length (%d) != Quality string length (%d)",
                     dna.id.c_str(), (int) dna.seq.size(), (int) dna.quals.size());
  for (size_t i = 0; i < dna.seq.size(); i++) {
    if ((dna.quals[i] - '!') < minimum_quality_score)
      dna.seq[i] = 'x';
  }
}

void ParseCommandLine(int argc, char **argv, Options &opts) {
  int opt;

  while ((opt = getopt(argc, argv, "h?H:t:o:T:p:R:C:U:O:Q:nmzqPSM")) != -1) {
    switch (opt) {
      case 'h' : case '?' :
        usage(0);
        break;
      case 'H' :
        opts.index_filename = optarg;
        break;
      case 't' :
        opts.taxonomy_filename = optarg;
        break;
      case 'T' :
        opts.confidence_threshold = std::stod(optarg);
        if (opts.confidence_threshold < 0 || opts.confidence_threshold > 1) {
          errx(EX_USAGE, "confidence threshold must be in [0, 1]");
        }
        break;
      case 'o' :
        opts.options_filename = optarg;
        break;
      case 'q' :
        opts.quick_mode = true;
        break;
      case 'p' :
        opts.num_threads = atoi(optarg);
        if (opts.num_threads < 1)
          errx(EX_USAGE, "number of threads can't be less than 1");
        break;
      case 'P' :
        opts.paired_end_processing = true;
        break;
      case 'S' :
        opts.paired_end_processing = true;
        opts.single_file_pairs = true;
        break;
      case 'm' :
        opts.mpa_style_report = true;
        break;
      case 'R' :
        opts.report_filename = optarg;
        break;
      case 'z' :
        opts.report_zero_counts = true;
        break;
      case 'C' :
        opts.classified_output_prefix = optarg;
        break;
      case 'U' :
        opts.unclassified_output_prefix = optarg;
        break;
      case 'O' :
        opts.kraken_output_filename = optarg;
        break;
      case 'n' :
        opts.print_scientific_name = true;
        break;
      case 'Q' :
        opts.minimum_quality_score = atoi(optarg);
        break;
      case 'M' :
        opts.use_memory_mapping = true;
        break;
    }
  }

  if (opts.index_filename.empty() ||
      opts.taxonomy_filename.empty() ||
      opts.options_filename.empty())
  {
    warnx("mandatory filename missing");
    usage();
  }

  if (opts.mpa_style_report && opts.report_filename.empty()) {
    warnx("-m requires -R be used");
    usage();
  }
}

void usage(int exit_code) {
  cerr << "Usage: classify [options] <fasta/fastq file(s)>" << endl
       << endl
       << "Options: (*mandatory)" << endl
       << "* -H filename      Kraken 2 index filename" << endl
       << "* -t filename      Kraken 2 taxonomy filename" << endl
       << "* -o filename      Kraken 2 options filename" << endl
       << "  -q               Quick mode" << endl
       << "  -M               Use memory mapping to access hash & taxonomy" << endl
       << "  -T NUM           Confidence score threshold (def. 0)" << endl
       << "  -p NUM           Number of threads (def. 1)" << endl
       << "  -Q NUM           Minimum quality score (FASTQ only, def. 0)" << endl
       << "  -P               Process pairs of reads" << endl
       << "  -S               Process pairs with mates in same file" << endl
       << "  -R filename      Print report to filename" << endl
       << "  -m               In comb. w/ -R, use mpa-style report" << endl
       << "  -z               In comb. w/ -R, report taxa w/ 0 count" << endl
       << "  -n               Print scientific name instead of taxid in Kraken output" << endl
       << "  -C prefix        Prefix for file with classified sequences" << endl
       << "  -U prefix        Prefix for file with unclassified sequences" << endl
       << "  -O filename      Output file for normal Kraken output" << endl;
  exit(exit_code);
}
