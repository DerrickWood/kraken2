/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include <err.h>
#include <sys/wait.h>
#include <sys/types.h>

#include "kraken2_headers.h"
#include "kv_store.h"
#include "taxonomy.h"
#include "seqreader.h"
#include "fast_reader.h"
#include "mmscanner.h"
#include "compact_hash.h"
#include "kraken2_data.h"
#include "aa_translate.h"
#include "reports.h"
#include "utilities.h"
using namespace kraken2;

using std::cerr;
using std::endl;
using std::ifstream;
using std::map;
using std::ostringstream;
using std::ofstream;
using std::string;
using std::vector;
using namespace kraken2;

static const size_t NUM_FRAGMENTS_PER_THREAD = 10000;
static const size_t INPUT_BLOCK_BYTES = 8 * 1024 * 1024;

// Mate identifiers agree once any trailing /1 or /2 is discounted.
static bool MatesAgree(const SeqView &a, const SeqView &b) {
  uint32_t la = a.header_len, lb = b.header_len;
  if (la > 2 && a.header[la - 2] == '/' &&
      (a.header[la - 1] == '1' || a.header[la - 1] == '2')) la -= 2;
  if (lb > 2 && b.header[lb - 2] == '/' &&
      (b.header[lb - 1] == '1' || b.header[lb - 1] == '2')) lb -= 2;
  return la == lb && memcmp(a.header, b.header, la) == 0;
}

// Re-emits a record from its view, appending a suffix to the header line.
static void WriteSeqView(ostringstream &oss, const SeqView &v,
                         const char *header_suffix) {
  oss << (v.format == FORMAT_FASTQ ? '@' : '>');
  oss.write(v.header, v.header_len);
  if (v.comment_len) {
    oss << ' ';
    oss.write(v.comment, v.comment_len);
  }
  oss << header_suffix << "\n";
  oss.write(v.seq, v.seq_len);
  oss << "\n";
  if (v.format == FORMAT_FASTQ) {
    oss << "+\n";
    oss.write(v.quals, v.quals_len);
    oss << "\n";
  }
}

// Masks in place over the reader's buffer.
static void MaskLowQualityBases(const SeqView &v, int minimum_quality_score) {
  if (v.quals == nullptr)
    return;
  char *seq = v.seq_mutable();
  for (uint32_t i = 0; i < v.seq_len; i++)
    if ((v.quals[i] - '!') < minimum_quality_score)
      seq[i] = 'x';
}
static const taxid_t MATE_PAIR_BORDER_TAXON = TAXID_MAX;
static const taxid_t READING_FRAME_BORDER_TAXON = TAXID_MAX - 1;
static const taxid_t AMBIGUOUS_SPAN_TAXON = TAXID_MAX - 2;

// Token stream used to defer minimizer lookups in ClassifySequence so they can be
// resolved in one prefetched batch. kind selects how replay rebuilds taxa[]/counts;
// key_idx indexes the distinct-minimizer lookup arrays for TOK_LOOKUP. uint32_t is
// ample: it counts distinct minimizers in a single read pair (far below 2^32).
enum MinTokKind { TOK_LOOKUP, TOK_SKIP, TOK_REPEAT, TOK_AMBIG,
                  TOK_BORDER_MATE, TOK_BORDER_FRAME };
struct MinToken { uint8_t kind; uint32_t key_idx; };

// template <typename Cell>
// using IndexData = std::tuple<IndexOptions, Taxonomy, CompactHashTable<Cell> *>;

struct IndexData {
  IndexOptions options;
  Taxonomy *taxonomy;
  KeyValueStore *cht;

  ~IndexData() {
    if (taxonomy) {
      delete taxonomy;
    }

    if (cht) {
      delete cht;
    }
  }
  // size_t cell_size;
};

struct Options {
  string index_filename;
  string taxonomy_filename;
  string options_filename;
  string report_filename;
  string classified_output_filename;
  string unclassified_output_filename;
  string kraken_output_filename;
  string taxon_counters_dump_filename;
  bool mpa_style_report;
  bool report_kmer_data;
  bool quick_mode;
  bool report_zero_counts;
  bool use_translated_search;
  bool print_scientific_name;
  double confidence_threshold;
  int num_threads;
  bool paired_end_processing;
  bool single_file_pairs;
  int minimum_quality_score;
  int minimum_hit_groups;
  bool use_memory_mapping;
  bool match_input_order;
  std::vector<char *> filenames;
  bool daemon_mode;
  bool check_pair_order;

  void reset() {
    quick_mode = false;
    confidence_threshold = 0;
    paired_end_processing = false;
    single_file_pairs = false;
    num_threads = 1;
    mpa_style_report = false;
    report_kmer_data = false;
    report_zero_counts = false;
    use_translated_search = false;
    print_scientific_name = false;
    minimum_quality_score = 0;
    minimum_hit_groups = 0;
    use_memory_mapping = false;
    daemon_mode = false;
    check_pair_order = false;

    index_filename.clear();
    taxonomy_filename.clear();
    options_filename.clear();
    report_filename.clear();
    classified_output_filename.clear();
    unclassified_output_filename.clear();
    kraken_output_filename.clear();
    filenames.clear();
    taxon_counters_dump_filename.clear();
  }
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
public:
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
    OutputStreamData &outputs, taxon_counters_t &total_taxon_counters);
taxid_t ClassifySequence(const SeqView &dna, const SeqView &dna2, ostringstream &koss,
    KeyValueStore *hash, Taxonomy &tax, IndexOptions &idx_opts,
    Options &opts, ClassificationStats &stats, MinimizerScanner &scanner,
    vector<taxid_t> &taxa, taxon_counts_t &hit_counts,
    vector<string> &tx_frames, taxon_counters_t &my_taxon_counts);
void AddHitlistString(ostringstream &oss, vector<taxid_t> &taxa,
    Taxonomy &taxonomy);
taxid_t ResolveTree(taxon_counts_t &hit_counts,
    Taxonomy &tax, size_t total_minimizers, Options &opts);
void ReportStats(struct timeval time1, struct timeval time2,
    ClassificationStats &stats);
void InitializeOutputs(Options &opts, OutputStreamData &outputs, SequenceFormat format);



void RemoveBlocking(int fd) {
  int flags = fcntl(fd, F_GETFL, 0);
  fcntl(fd, F_SETFL, flags & ~O_NONBLOCK);
}

int daemonize() {
  pid_t pid;

  if ((pid = fork()) < 0) {
    errx(1, "fork");
  }

  if (pid == 0) {
    if (setsid() == -1) {
      perror("setdsid");
      exit(1);
    }
  } else {
    exit(0);
  }

  if ((pid = fork()) < 0) {
    perror("fork");
    exit(1);
  }

  if (pid != 0) {
    exit(0);
  }

  mkfifo("/tmp/classify_stdin", S_IRWXU);
  mkfifo("/tmp/classify_stdout", S_IRWXU);

  int read_fd = open("/tmp/classify_stdin", O_RDONLY | O_NONBLOCK);
  // keep these open so that to daemon does not
  // receive EOF when the external process closes
  // its end of the fifo
  int dummy_fd_1 = open("/tmp/classify_stdin", O_WRONLY);
  (void)dummy_fd_1;
  int dummy_fd_2 = open("/tmp/classify_stdout", O_RDONLY | O_NONBLOCK);

  int write_fd = open("/tmp/classify_stdout", O_WRONLY);

  RemoveBlocking(read_fd);
  RemoveBlocking(dummy_fd_2);

  for (int fd = 0; fd < 2; fd++) {
    close(fd);
  }

  dup2(read_fd, 0);
  dup2(write_fd, 1);
  dup2(write_fd, 2);

  return (0);
}

void OpenFifos(Options &opts, pid_t pid) {
  char stdin_filename[256];
  char stdout_filename[256];

  snprintf(stdin_filename, 256, "/tmp/classify_%d_stdin", pid);
  snprintf(stdout_filename, 256, "/tmp/classify_%d_stdout", pid);

  mkfifo(stdin_filename, S_IRWXU);
  mkfifo(stdout_filename, S_IRWXU);

  int read_fd = -1;
  // if we are expecting input from stdin, open the
  // FIFO in blocking mode and wait for the input.
  if (opts.filenames.empty()) {
    read_fd = open(stdin_filename, O_RDONLY);
  } else {
    read_fd = open(stdin_filename, O_RDONLY | O_NONBLOCK);
    // keep these open so that to daemon does not
    // receive EOF when the external process closes
    // its end of the fifo
    int dummy_fd_1 = open(stdin_filename, O_WRONLY);
    (void)dummy_fd_1;
  }

  // If we are outputting to a file open the write
  // FIFO in non-blocking mode and keep a read end open
  // so that we do not block the process.
  int dummy_fd_2 = -1;
  if (!opts.kraken_output_filename.empty()) {
    dummy_fd_2 = open(stdout_filename, O_RDONLY | O_NONBLOCK);
  }
  int write_fd = open(stdout_filename, O_WRONLY);

  if (opts.kraken_output_filename.empty()) {
    RemoveBlocking(read_fd);
  } else {
    RemoveBlocking(dummy_fd_2);
  }

  for (int fd = 0; fd < 3; fd++) {
    close(fd);
  }

  dup2(read_fd, 0);
  dup2(write_fd, 1);
  dup2(write_fd, 2);
}

IndexData *load_index(Options &opts) {
  cerr << "Loading database information...";

  IndexOptions idx_opts = {0};
  ifstream idx_opt_fs(opts.options_filename);
  struct stat sb;
  if (stat(opts.options_filename.c_str(), &sb) < 0)
    errx(EX_OSERR, "unable to get filesize of %s", opts.options_filename.c_str());
  auto opts_filesize = sb.st_size;
  idx_opt_fs.read((char *) &idx_opts, opts_filesize);
  opts.use_translated_search = ! idx_opts.dna_db;

  IndexData *index_data = new IndexData();
  index_data->options = idx_opts;
  index_data->taxonomy = new Taxonomy(opts.taxonomy_filename, opts.use_memory_mapping);

  KeyValueStore *cht;
  switch (GetKVStoreCellType(opts.index_filename)) {
  case CompactHash32:
    cht = new CompactHashTable<CompactHashCell>(opts.index_filename, opts.use_memory_mapping);
    break;
  case CompactHash40:
    cht = new CompactHashTable<CompactHashCell40>(opts.index_filename, opts.use_memory_mapping);
    break;
  case Unknown:
    errx(1, "Unable to determine width of compact hash cell\n");
  }

  index_data->cht = cht;

  cerr << " done." << endl;

  return index_data;
}

void classify(Options &opts, IndexData *index_data) {
  taxon_counters_t taxon_counters; // stats per taxon
  IndexOptions idx_opts = index_data->options;
  Taxonomy &taxonomy = *(index_data->taxonomy);
  KeyValueStore *hash_ptr = index_data->cht;

  omp_set_num_threads(opts.num_threads);

  ClassificationStats stats = {0, 0, 0};

  OutputStreamData outputs = { false, false, nullptr, nullptr, nullptr, nullptr, &std::cout };
  // taxon_counters.reserve(taxonomy.node_count());
  struct timeval tv1, tv2;
  gettimeofday(&tv1, nullptr);
  if (opts.filenames.empty()) {
    if (opts.paired_end_processing && ! opts.single_file_pairs)
      errx(EX_USAGE, "paired end processing used with no files specified");
    ProcessFiles(nullptr, nullptr, hash_ptr, taxonomy, idx_opts, opts, stats, outputs, taxon_counters);
  }
  else {
    for (size_t i = 0; i < opts.filenames.size(); i++) {
      if (opts.paired_end_processing && ! opts.single_file_pairs) {
        if (i + 1 == opts.filenames.size()) {
          errx(EX_USAGE, "paired end processing used with unpaired file");
        }
        ProcessFiles(opts.filenames[i], opts.filenames[i+1], hash_ptr, taxonomy, idx_opts, opts, stats, outputs, taxon_counters);
        i += 1;
      } else {
        ProcessFiles(opts.filenames[i], nullptr, hash_ptr, taxonomy, idx_opts, opts, stats, outputs, taxon_counters);
      }
    }
  }
  gettimeofday(&tv2, nullptr);

  // delete hash_ptr;

  ReportStats(tv1, tv2, stats);
  if (!opts.taxon_counters_dump_filename.empty()) {
    dumpTaxonCounters(opts.taxon_counters_dump_filename, taxon_counters, taxonomy);
  }
  if (! opts.report_filename.empty()) {
    if (opts.mpa_style_report)
      ReportMpaStyle(opts.report_filename, opts.report_zero_counts, taxonomy,
          taxon_counters);
    else {
      auto total_unclassified = stats.total_sequences - stats.total_classified;
      ReportKrakenStyle(opts.report_filename, opts.report_zero_counts,
          opts.report_kmer_data, taxonomy,
          taxon_counters, stats.total_sequences, total_unclassified);
    }
  }
}

void TokenizeString(std::vector<char *>& argv, const char *str) {
  argv.resize(0);
  const char *sep = " ";
  char *str_cpy = new char[strlen(str) + 1];
  char *token;

  strcpy(str_cpy, str);
  for (token = strtok(str_cpy, sep); token;
       token = strtok(NULL, sep)) {
    int token_len = strlen(token);
    char *word = NULL;
    word = new char[token_len + 1];
    strcpy(word, token);
    argv.push_back(word);
  }

  delete[] str_cpy;
}

void ClassifyDaemon(Options opts) {
  daemonize();

  std::vector<char *> args;
  std::map<std::string, IndexData*> indexes;
  std::stringstream ss;
  char *cmdline = NULL;
  size_t linecap = 0;
  ssize_t line_len;
  bool stop = false;

  int pid_file = open("/tmp/classify.pid", O_CREAT | O_WRONLY | O_TRUNC, 0600);
  ss << getpid() << "\n";
  write(pid_file, ss.str().c_str(), ss.str().size());
  ss.str("");
  close(pid_file);

  IndexData *index_data = load_index(opts);
  indexes[opts.index_filename] = index_data;

  while (!stop) {
    classify(opts, index_data);
    std::cout << "DONE" << std::endl;

    while (true) {
      line_len = getline(&cmdline, &linecap, stdin);
      if (line_len < 2)
        continue;
      if (cmdline == std::string("PING\n")) {
        std::cerr << "OK" << std::endl;
        continue;
      } else if (cmdline == std::string("STOP\n")) {
        std::cerr << "OK" << std::endl;
        stop = true;
      }
      cmdline[line_len - 1] = ' ';
      break;
    }

    TokenizeString(args, cmdline);
    opts.reset();
    ParseCommandLine(args.size(), &args[0], opts);
    if (indexes.find(opts.index_filename) == indexes.end()) {
      index_data = load_index(opts);
      indexes[opts.index_filename] = index_data;
    }
  }

  for (auto i = 0; i < 3; i++) {
    close(i);
  }

  free(cmdline);

  for (auto kv : indexes) {
    delete kv.second;  // ~IndexData() frees cht + taxonomy (was double-free)
  }
}

int main(int argc, char **argv) {
  Options opts;
  opts.reset();
  ParseCommandLine(argc, argv, opts);
  if (opts.daemon_mode) {
    ClassifyDaemon(opts);
  } else {
    IndexData *index_data = load_index(opts);
    classify(opts, index_data);
    delete index_data;  // ~IndexData() frees cht + taxonomy (was double-free)
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

  if (isatty(fileno(stderr)))
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
    OutputStreamData &outputs,
    taxon_counters_t &total_taxon_counters)
{
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
  // The critical section reads raw bytes from these descriptors and cuts on a
  // record boundary; parsing happens afterwards, outside the lock.
  int fd1 = filename1 ? open(filename1, O_RDONLY) : fileno(stdin);
  int fd2 = filename2 ? open(filename2, O_RDONLY) : -1;
  if (fd1 < 0)
    errx(EX_NOINPUT, "unable to open %s", filename1);
  if (filename2 && fd2 < 0)
    errx(EX_NOINPUT, "unable to open %s", filename2);
  StreamCursor cursor1, cursor2;
  std::vector<OutputData> buffers(omp_get_max_threads() * 2);

  #pragma omp parallel
  {
    MinimizerScanner scanner(idx_opts.k, idx_opts.l, idx_opts.spaced_seed_mask,
                             idx_opts.dna_db, idx_opts.toggle_mask,
                             idx_opts.revcom_version);
    vector<taxid_t> taxa;
    taxon_counts_t hit_counts;
    ostringstream kraken_oss, c1_oss, c2_oss, u1_oss, u2_oss;
    ClassificationStats thread_stats = {0, 0, 0};
    vector<string> translated_frames(6);
    SeqView *seq1 = nullptr, *seq2 = nullptr;
    FastReader reader1, reader2;
    size_t idx1 = 0, idx2 = 0;
    uint64_t block_id;
    OutputData out_data;
    taxon_counters_t thread_taxon_counters;

    while (true) {
      thread_stats.total_sequences = 0;
      thread_stats.total_bases = 0;
      thread_stats.total_classified = 0;

      auto ok_read = false;

      #pragma omp critical(seqread)
      {  // Input processing block
        if (! opts.paired_end_processing) {
          ok_read = reader1.LoadBlock(fd1, cursor1, INPUT_BLOCK_BYTES);
        }
        else if (! opts.single_file_pairs) {
          // Take a block from the first mate, then exactly as many records from
          // the second, so the two files stay in step.
          ok_read = reader1.LoadBlock(fd1, cursor1, INPUT_BLOCK_BYTES);
          if (ok_read)
            ok_read = reader2.LoadRecords(fd2, cursor2, reader1.RecordCount());
        }
        else {
          // Interleaved pairs: cut on an even record count so a pair is never
          // split across blocks.
          ok_read = reader1.LoadBlock(fd1, cursor1, INPUT_BLOCK_BYTES, 2);
        }
        block_id = next_input_block_id++;
      }

      if (! ok_read)
        break;

      // Parsing is deliberately outside the critical section above.
      reader1.Parse();
      if (opts.paired_end_processing && ! opts.single_file_pairs)
        reader2.Parse();
      idx1 = idx2 = 0;

      // printing_sequences gates whether records are emitted below, so the
      // outputs must be opened before the first block is classified.
      if (! outputs.initialized)
        InitializeOutputs(opts, outputs, reader1.file_format());

      // Reset all dynamically-growing things
      kraken_oss.str("");
      c1_oss.str("");
      c2_oss.str("");
      u1_oss.str("");
      u2_oss.str("");
      thread_taxon_counters.clear();

      while (idx1 < reader1.size()) {
        seq1 = &reader1.at(idx1++);
        auto valid_fragment = true;
        if (opts.paired_end_processing && valid_fragment) {
          if (opts.single_file_pairs) {
            valid_fragment = idx1 < reader1.size();
            seq2 = valid_fragment ? &reader1.at(idx1++) : nullptr;
          } else {
            valid_fragment = idx2 < reader2.size();
            seq2 = valid_fragment ? &reader2.at(idx2++) : nullptr;
          }
          if (! valid_fragment)
            break;
          if (! MatesAgree(*seq1, *seq2)) {
            errx(1, "ERROR: Unmatched pairs.\n"
                 "Mate 1: %.*s\nMate 2: %.*s.\nPlease make sure that pairs "
                 "are sorted before classification.",
                 (int) seq1->header_len, seq1->header,
                 (int) seq2->header_len, seq2->header);
          }
        }
        if (! valid_fragment)
          break;
        thread_stats.total_sequences++;
        if (opts.minimum_quality_score > 0) {
          MaskLowQualityBases(*seq1, opts.minimum_quality_score);
          if (opts.paired_end_processing)
            MaskLowQualityBases(*seq2, opts.minimum_quality_score);
        }
        taxid_t call;
        if (opts.paired_end_processing) {
          call =
              ClassifySequence(*seq1, *seq2, kraken_oss, hash, tax, idx_opts,
                               opts, thread_stats, scanner, taxa, hit_counts,
                               translated_frames, thread_taxon_counters);
        } else {
          static const SeqView empty_sequence = { nullptr, nullptr, nullptr,
              nullptr, 0, 0, 0, 0, FORMAT_FASTQ };
          call = ClassifySequence(*seq1, empty_sequence, kraken_oss, hash, tax, idx_opts,
                                  opts, thread_stats, scanner, taxa, hit_counts,
                                  translated_frames, thread_taxon_counters);
        }
        if (outputs.printing_sequences) {
          char buffer[64] = "";
          if (call)
            sprintf(buffer, " kraken:taxid|%llu",
                (unsigned long long) tax.nodes()[call].external_id);
          WriteSeqView(call ? c1_oss : u1_oss, *seq1, call ? buffer : "");
          if (opts.paired_end_processing)
            WriteSeqView(call ? c2_oss : u2_oss, *seq2, call ? buffer : "");
        }
        thread_stats.total_bases += seq1->seq_len;
        if (opts.paired_end_processing)
          thread_stats.total_bases += seq2->seq_len;
      }

      // #pragma omp atomic
      // #pragma omp atomic
      // #pragma omp atomic

      #pragma omp critical(output_stats)
      {
        stats.total_bases += thread_stats.total_bases;
        stats.total_sequences += thread_stats.total_sequences;
        stats.total_classified += thread_stats.total_classified;

        if (isatty(fileno(stderr)))
          cerr << "\rProcessed " << stats.total_sequences
               << " sequences (" << stats.total_bases << " bp) ...";
      }

      out_data.block_id = block_id;
      out_data.kraken_str.assign(kraken_oss.str());
      out_data.classified_out1_str.assign(c1_oss.str());
      out_data.classified_out2_str.assign(c2_oss.str());
      out_data.unclassified_out1_str.assign(u1_oss.str());
      out_data.unclassified_out2_str.assign(u2_oss.str());

      #pragma omp critical(output_queue)
      {
        output_queue.push(std::move(out_data));
      }

      if (!opts.report_filename.empty() || !opts.taxon_counters_dump_filename.empty()) {
#pragma omp critical(update_taxon_counters)
        for (auto &kv_pair : thread_taxon_counters) {
          total_taxon_counters[kv_pair.first] += std::move(kv_pair.second);
        }
      }

      bool output_loop = true;
      bool borrowed_buffer = false;
      while (output_loop) {
        #pragma omp critical(output_queue)
        {
          output_loop = !output_queue.empty();
          if (!output_loop && !borrowed_buffer) {
            out_data = std::move(buffers.back());
            buffers.pop_back();
          } else if (borrowed_buffer && output_loop) {
            buffers.push_back(std::move(out_data));
          }

          if (output_loop) {
            if (output_queue.top().block_id == next_output_block_id) {
              out_data = std::move(const_cast<OutputData &>(output_queue.top()));
              output_queue.pop();
              borrowed_buffer = true;
              // Acquiring output lock obligates thread to print out
              // next output data block, contained in out_data
              omp_set_lock(&output_lock);
              next_output_block_id++;
            } else {
              output_loop = false;
              if (buffers.size() > 0) {
                out_data = std::move(buffers.back());
                buffers.pop_back();
              } else {
                out_data = OutputData();
              }
            }
          }
        }
        if (! output_loop)
          break;
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
    } // end while
  } // end parallel block

  omp_destroy_lock(&output_lock);
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

  // Sum each taxon's LTR path, find taxon with highest LTR score
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

  // Reset max. score to be only hits at the called taxon
  max_score = hit_counts[max_taxon];
  // We probably have a call w/o required support (unless LCA resolved tie)
  while (max_taxon && max_score < required_score) {
    max_score = 0;
    for (auto &kv_pair : hit_counts) {
      taxid_t taxon = kv_pair.first;
      // Add to score if taxon in max_taxon's clade
      if (taxonomy.IsAAncestorOfB(max_taxon, taxon)) {
        max_score += kv_pair.second;
      }
    }
    // Score is now sum of hits at max_taxon and w/in max_taxon clade
    if (max_score >= required_score)
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
  if ( id[sz - 2] == '/' && (id[sz - 1] == '1' || id[sz - 1] == '2') )
    return id.substr(0, sz - 2);
  return id;
}

taxid_t ClassifySequence(const SeqView &dna, const SeqView &dna2, ostringstream &koss,
                         KeyValueStore *hash, Taxonomy &taxonomy,
                         IndexOptions &idx_opts, Options &opts,
                         ClassificationStats &stats, MinimizerScanner &scanner,
                         vector<taxid_t> &taxa, taxon_counts_t &hit_counts,
                         vector<string> &tx_frames,
                         taxon_counters_t &curr_taxon_counts)
{
  uint64_t *minimizer_ptr;
  taxid_t call = 0;
  taxa.clear();
  hit_counts.clear();
  auto frame_ct = opts.use_translated_search ? 6 : 1;
  int64_t minimizer_hit_groups = 0;

  // Phase 1: scan minimizers for the read pair into a token stream, deferring the
  // memory-latency-bound hash lookups. Lookup-vs-repeat-vs-ambiguous is decided
  // purely from minimizer values (independent of taxon), exactly as the original
  // consecutive-dedup did, so phase 3 reproduces identical taxa[]/counts.
  static thread_local std::vector<uint64_t> lookup_keys;
  static thread_local std::vector<MinToken> tok_stream;
  lookup_keys.clear();
  tok_stream.clear();

  for (int mate_num = 0; mate_num < 2; mate_num++) {
    if (mate_num == 1 && ! opts.paired_end_processing)
      break;

    const SeqView &mate = (mate_num == 0) ? dna : dna2;
    if (opts.use_translated_search) {
      // The frame translator wants a real string; this path is the rare one.
      std::string mate_seq(mate.seq, mate.seq_len);
      TranslateToAllFrames(mate_seq, tx_frames);
    }
    // index of frame is 0 - 5 w/ tx search (or 0 if no tx search)
    for (int frame_idx = 0; frame_idx < frame_ct; frame_idx++) {
      if (opts.use_translated_search) {
        scanner.LoadSequence(tx_frames[frame_idx]);
      }
      else {
        scanner.LoadSequence(mate.seq, mate.seq_len);
      }
      uint64_t last_minimizer = UINT64_MAX;
      while ((minimizer_ptr = scanner.NextMinimizer()) != nullptr) {
        if (scanner.is_ambiguous()) {
          tok_stream.push_back({TOK_AMBIG, 0});
        }
        else if (*minimizer_ptr != last_minimizer) {
          last_minimizer = *minimizer_ptr;
          bool skip_lookup = idx_opts.minimum_acceptable_hash_value &&
              MurmurHash3(*minimizer_ptr) < idx_opts.minimum_acceptable_hash_value;
          if (skip_lookup) {
            tok_stream.push_back({TOK_SKIP, 0});
          }
          else {
            tok_stream.push_back({TOK_LOOKUP, (uint32_t) lookup_keys.size()});
            lookup_keys.push_back(*minimizer_ptr);
          }
        }
        else {
          tok_stream.push_back({TOK_REPEAT, 0});
        }
      }
      if (opts.use_translated_search && frame_idx != 5)
        tok_stream.push_back({TOK_BORDER_FRAME, 0});
    }
    if (opts.paired_end_processing && mate_num == 0)
      tok_stream.push_back({TOK_BORDER_MATE, 0});
  }

  // Phase 2: resolve all distinct minimizers in one prefetched batched pass.
  static thread_local std::vector<hvalue_t> lookup_vals;
  lookup_vals.resize(lookup_keys.size());
  if (! lookup_keys.empty())
    hash->GetBatch(lookup_keys.data(), lookup_vals.data(), lookup_keys.size());

  // Phase 3: replay token stream, reproducing the original behavior exactly.
  {
    taxid_t last_taxon = 0;
    for (size_t ti = 0; ti < tok_stream.size(); ti++) {
      const MinToken &tok = tok_stream[ti];
      taxid_t taxon = 0;
      switch (tok.kind) {
        case TOK_AMBIG:
          taxa.push_back(AMBIGUOUS_SPAN_TAXON);
          continue;
        case TOK_BORDER_FRAME:
          taxa.push_back(READING_FRAME_BORDER_TAXON);
          continue;
        case TOK_BORDER_MATE:
          taxa.push_back(MATE_PAIR_BORDER_TAXON);
          continue;
        case TOK_SKIP:
          taxon = 0;
          last_taxon = 0;
          break;
        case TOK_LOOKUP:
          taxon = lookup_vals[tok.key_idx];
          last_taxon = taxon;
          if (taxon) {
            minimizer_hit_groups++;
            if (!opts.report_filename.empty() || !opts.taxon_counters_dump_filename.empty())
              curr_taxon_counts[taxon].add_kmer(lookup_keys[tok.key_idx]);
          }
          break;
        default:  // TOK_REPEAT
          taxon = last_taxon;
          curr_taxon_counts[taxon].add_kmer(lookup_keys[tok.key_idx]);
          break;
      }
      if (taxon) {
        if (opts.quick_mode && minimizer_hit_groups >= opts.minimum_hit_groups) {
          call = taxon;
          goto finished_searching;
        }
        hit_counts[taxon]++;
      }
      taxa.push_back(taxon);
    }
  }

  finished_searching:

  auto total_kmers = taxa.size();
  if (opts.paired_end_processing)
    total_kmers--;  // account for the mate pair marker
  if (opts.use_translated_search)  // account for reading frame markers
    total_kmers -= opts.paired_end_processing ? 4 : 2;
  call = ResolveTree(hit_counts, taxonomy, total_kmers, opts);
  // Void a call made by too few minimizer groups
  if (call && minimizer_hit_groups < opts.minimum_hit_groups)
    call = 0;

  if (call) {
    stats.total_classified++;
    if (!opts.report_filename.empty() || !opts.taxon_counters_dump_filename.empty())
      curr_taxon_counts[call].incrementReadCount();
  }

  if (call)
    koss << "C\t";
  else
    koss << "U\t";
  {
    uint32_t n = dna.header_len;
    if (opts.paired_end_processing && n > 2 && dna.header[n - 2] == '/' &&
        (dna.header[n - 1] == '1' || dna.header[n - 1] == '2'))
      n -= 2;
    koss.write(dna.header, n);
    koss << "\t";
  }

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
    koss << dna.seq_len << "\t";
  else
    koss << dna.seq_len << "|" << dna2.seq_len << "\t";

  if (opts.quick_mode) {
    koss << ext_call << ":Q";
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

ofstream *OpenOutputStream(const std::string &filename) {
  ofstream *out = new ofstream();
  out->exceptions(std::ios::badbit | std::ios::failbit);
  try {
    out->open(filename);
  } catch(std::exception &e) {
    std::cerr << '\r' << "Unable to open file: " << filename
              << ", reason: " << strerror(errno) << std::endl;
    delete out;
    exit(EXIT_FAILURE);
  }

  return out;
}

void InitializeOutputs(Options &opts, OutputStreamData &outputs, SequenceFormat format) {
  #pragma omp critical(output_init)
  {
    if (! outputs.initialized) {
      if (! opts.classified_output_filename.empty()) {
        if (opts.paired_end_processing) {
          vector<string> fields = SplitString(opts.classified_output_filename, "#", 3);
          if (fields.size() < 2) {
            errx(EX_DATAERR, "Paired filename format missing # character: %s",
                 opts.classified_output_filename.c_str());
          }
          else if (fields.size() > 2) {
            errx(EX_DATAERR, "Paired filename format has >1 # character: %s",
                 opts.classified_output_filename.c_str());
          }
          outputs.classified_output1 = OpenOutputStream(fields[0] + "_1" + fields[1]);
          outputs.classified_output2 = OpenOutputStream(fields[0] + "_2" + fields[1]);
        }
        else
          outputs.classified_output1 = new ofstream(opts.classified_output_filename);
        outputs.printing_sequences = true;
      }
      if (! opts.unclassified_output_filename.empty()) {
        if (opts.paired_end_processing) {
          vector<string> fields = SplitString(opts.unclassified_output_filename, "#", 3);
          if (fields.size() < 2) {
            errx(EX_DATAERR, "Paired filename format missing # character: %s",
                 opts.unclassified_output_filename.c_str());
          }
          else if (fields.size() > 2) {
            errx(EX_DATAERR, "Paired filename format has >1 # character: %s",
                 opts.unclassified_output_filename.c_str());
          }
          outputs.unclassified_output1 = OpenOutputStream(fields[0] + "_1" + fields[1]);
          outputs.unclassified_output2 = OpenOutputStream(fields[0] + "_2" + fields[1]);
        }
        else
          outputs.unclassified_output1 = new ofstream(opts.unclassified_output_filename);
        outputs.printing_sequences = true;
      }
      if (!opts.kraken_output_filename.empty()) {
        if (opts.kraken_output_filename == "-")  // Special filename to silence Kraken output
          outputs.kraken_output = nullptr;
        else {
          // outputs.kraken_output = new ofstream(opts.kraken_output_filename);
          outputs.kraken_output = OpenOutputStream(opts.kraken_output_filename);
        }
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
                     dna.header.c_str(), (int) dna.seq.size(), (int) dna.quals.size());
  for (size_t i = 0; i < dna.seq.size(); i++) {
    if ((dna.quals[i] - '!') < minimum_quality_score)
      dna.seq[i] = 'x';
  }
}

void ParseCommandLine(int argc, char **argv, Options &opts) {
  int opt;

  while ((opt = getopt(argc, argv, "h?H:t:o:T:p:R:C:U:O:Q:g:d:nmzqPSMKDc")) != -1) {
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
      case 'g' :
        opts.minimum_hit_groups = atoi(optarg);
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
      case 'K':
        opts.report_kmer_data = true;
        break;
      case 'R' :
        opts.report_filename = optarg;
        break;
      case 'z' :
        opts.report_zero_counts = true;
        break;
      case 'C' :
        opts.classified_output_filename = optarg;
        break;
      case 'U' :
        opts.unclassified_output_filename = optarg;
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
      case 'D':
        opts.daemon_mode = true;
        break;
      case 'c':
        opts.check_pair_order = true;
        break;
      case 'd':
        opts.taxon_counters_dump_filename = optarg;
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

  for (int i = optind; i < argc; i++) {
    opts.filenames.push_back(argv[i]);
  }

  optind = 1;
}

void usage(int exit_code) {
  cerr << "Usage: classify [options] <fasta/fastq file(s)>" << endl
       << endl
       << "Options: (*mandatory)" << endl
       << "* -H filename      Kraken 2 index filename" << endl
       << "* -t filename      Kraken 2 taxonomy filename" << endl
       << "* -o filename      Kraken 2 options filename" << endl
       << "  -q               Quick mode" << endl
       << "  -c               Ensure pairs are ordered (stop classification otherwise)"
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
       << "  -g NUM           Minimum number of hit groups needed for call" << endl
       << "  -C filename      Filename/format to have classified sequences" << endl
       << "  -U filename      Filename/format to have unclassified sequences" << endl
       << "  -O filename      Output file for normal Kraken output" << endl
       << "  -K               In comb. w/ -R, provide minimizer information in report" << endl
       << "  -D               Start a daemon, this options is intended to be used with wrappers" << std::endl
       << "  -d filename      Dump taxon counters to filename." << endl;
    exit(exit_code);
}
