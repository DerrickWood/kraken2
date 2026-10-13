/*
 * Symmetric Dustmasker mostly based on a similarly named implementation
 * by Heng Li as part of the minimap project.
 */

#include <algorithm>
#include <deque>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include <queue>

#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include "fast_reader.h"

using namespace kraken2;

// input_lock serializes reads from the input stream. output_lock guards the
// output stream, pending_blocks and next_output_block, so blocks are written
// in input order.
static std::mutex input_lock;
static std::mutex output_lock;
static std::priority_queue<size_t, std::vector<size_t>, std::greater<size_t>> pending_blocks;
static size_t next_output_block = 0;

// Directory for blocks that finish masking before the blocks ahead of them
// have been written. Spilling them to disk bounds memory use.
static std::string spill_dir;

uint8_t asc2dna[] = {
  /*   0 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  /*  16 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  /*  32 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  /*                                               - */
  /*  48 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  /*  64 */ 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
         /*    A  B  C  D        G  H        K     M  N */
  /*  80 */ 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
         /*       R  S  T  U  V  W     Y */
  /*  96 */ 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
         /*    a  b  c  d        g  h        k     m  n */
  /* 112 */ 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
         /*       r  s  t  u  v  w     y */
  /* 128 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  /* 144 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  /* 160 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  /* 176 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  /* 192 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  /* 208 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  /* 224 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  /* 240 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
};

static const size_t INPUT_BLOCK_BYTES = 8 * 1024 * 1024;

struct PerfectInterval {
  int start;
  int finish;
  int left;
  int right;
};

struct MaskRange {
  int start;
  int finish;
};

struct SDust {
  SDust() : rw(0), rv(0), l(0) {
    memset(cv, 0, 64 * sizeof(int));
    memset(cw, 0, 64 * sizeof(int));
  }

  void reset() {
    kmers.clear();
    perfectIntervals.clear();
    ranges.clear();
    memset(cw, 0, 64 * sizeof(int));
    memset(cv, 0, 64 * sizeof(int));
    l = rw = rv = 0;
  }

  std::deque<int> kmers;
  std::vector<PerfectInterval> perfectIntervals;
  std::vector<MaskRange> ranges;
  SeqView *seq;
  int cw[64];
  int cv[64];
  int rw;
  int rv;
  int l;
};

int windowSize = 64;
int threshold = 20;
std::function<int(int)> processMaskedNucleotide = tolower;

void usage(const char *prog) {
  std::cerr << "usage: " << prog
            << " [-T | -level threshold] [-W | -window window size] [-i | -in "
               "input file]"
            << std::endl
            << " [-w | -width sequence width] [-o | -out output file] [-r | "
               "-replace-masked-with char]"
            << std::endl
            << " [-f | -outfmt output format] [-t | -threads threads]"
            << std::endl;
  exit(1);
}

bool fileExists(const char *filename) {
  std::ifstream ifs(filename);

  return ifs.good();
}

bool stricasecmp(std::string &s1, std::string &s2) {
  if (s1.size() != s2.size())
    return false;
  for (size_t i = 0; i < s1.size(); i++) {
    if (tolower(s1[i]) != tolower(s2[i]))
      return false;
  }
  return true;
}

void shiftWindow(SDust &sd, int t) {
  int s;
  if ((int)sd.kmers.size() >= windowSize - 2) {
    s = sd.kmers.front();
    sd.kmers.pop_front();
    sd.rw -= --sd.cw[s];
    if (sd.l > (int)sd.kmers.size()) {
      --sd.l;
      sd.rv -= --sd.cv[s];
    }
  }
  sd.kmers.push_back(t);
  sd.l++;
  sd.rw += sd.cw[t]++;
  sd.rv += sd.cv[t]++;
  if (sd.cv[t] * 10 > threshold * 2) {
    do {
      s = sd.kmers.at(sd.kmers.size() - sd.l);
      sd.rv -= --sd.cv[s];
      sd.l--;
    } while (s != t);
  }
}

void saveMaskedRegions(SDust &sd, int windowStart) {
  bool saved = false;
  if (sd.perfectIntervals.size() == 0 ||
      sd.perfectIntervals.back().start >= windowStart)
    return;
  PerfectInterval &p = sd.perfectIntervals.back();
  if (!sd.ranges.empty()) {
    int start = sd.ranges.back().start;
    int finish = sd.ranges.back().finish;
    if (p.start <= finish) {
      sd.ranges.back() = MaskRange{start, std::max(p.finish, finish)};
      saved = true;
    }
  }
  if (!saved)
    sd.ranges.push_back(MaskRange{p.start, p.finish});
  while (!sd.perfectIntervals.empty() &&
         sd.perfectIntervals.back().start < windowStart)
    sd.perfectIntervals.pop_back();
}

void findPerfect(SDust &sd, int windowStart) {
  int cv[64];
  int maxLeft = 0;
  int maxRight = 0;
  int newLeft = 0;
  int newRight = sd.rv;

  memcpy(cv, sd.cv, 64 * sizeof(int));
  for (int i = sd.kmers.size() - sd.l - 1; i >= 0; --i) {
    size_t j;
    int kmer = sd.kmers.at(i);
    newRight += cv[kmer]++;
    newLeft = sd.kmers.size() - i - 1;
    if (newRight * 10 > threshold * newLeft) {
      for (j = 0; j < sd.perfectIntervals.size() &&
                  sd.perfectIntervals[j].start >= i + windowStart;
           ++j) {
        PerfectInterval &p = sd.perfectIntervals[j];
        if (maxRight == 0 || p.right * maxLeft > maxRight * p.left) {
          maxLeft = p.left;
          maxRight = p.right;
        }
      }
      if (maxRight == 0 || newRight * maxLeft >= maxRight * newLeft) {
        maxLeft = newLeft;
        maxRight = newRight;
        PerfectInterval p;
        p.start = i + windowStart;
        p.finish = sd.kmers.size() + 2 + windowStart;
        p.left = newLeft;
        p.right = newRight;
        auto position = sd.perfectIntervals.begin() + j;
        sd.perfectIntervals.insert(position, p);
      }
    }
  }
}

void runSymmetricDust(SDust &sd, char *seq, size_t size, int offset) {
  int triplet = 0;
  int windowStart = 0;
  int l = 0;
  for (size_t i = 0; i < size; i++) {
    int base = asc2dna[(int)seq[i]];
    if (base < 4) {
      l++;
      triplet = (triplet << 2 | base) & 63;
      if (l >= 3) {
        windowStart = std::max(l - (int)windowSize, 0);
        saveMaskedRegions(sd, windowStart);
        shiftWindow(sd, triplet);
        if (sd.rw * 10 > sd.l * threshold)
          findPerfect(sd, windowStart);
      }
    }
  }
  while (!sd.perfectIntervals.empty())
    saveMaskedRegions(sd, windowStart++);
}

void printSequence(SeqView *seq, std::ostream &out, int width) {
  out << ">";
  out.write(seq->header, seq->header_len);
  if (seq->comment_len) {
    out << " ";
    out.write(seq->comment, seq->comment_len);
  }
  out << "\n";
  for (size_t i = 0; i < seq->seq_len; i += width) {
    if (i + width >= seq->seq_len)
      width = seq->seq_len - i;
    out.write(seq->seq + i, width);
    out << '\n';
  }
}

void mask(SDust &sd) {
  size_t i = 0;
  char *seq = sd.seq->seq_mutable();
  size_t seq_len = sd.seq->seq_len;

  for (; i < seq_len; i++) {
    if (asc2dna[(int)seq[i]] != 4) {
      int start = i;
      for (;;) {
        seq[i] = toupper(seq[i]);
        if ((i + 1) == seq_len || asc2dna[(int)seq[i + 1]] == 4)
          break;
        i++;
      }
      char *s = &seq[0] + start;
      runSymmetricDust(sd, s, i - start + 1, start);
      for (size_t j = 0; j < sd.ranges.size(); j++) {
        for (int i = sd.ranges[j].start; i < sd.ranges[j].finish; i++)
          s[i] = processMaskedNucleotide(s[i]);
      }
      sd.reset();
    }
  }
}

void flushOutput(FastReader &reader, std::ostream &out, size_t width) {
  for (size_t i = 0; i < reader.size(); i++) {
    printSequence(&reader.at(i), out, width);
  }
}

std::string spillPath(size_t block_number) {
  return spill_dir + "/block." + std::to_string(block_number);
}

void fail(const std::string &what) {
  std::cerr << "k2mask: " << what << ": " << strerror(errno) << std::endl;
  std::exit(1);
}

// Writes every spilled block that is next in input order.  The caller holds
// output_lock.
void writeReadySpilledBlocks(std::ostream &out) {
  while (!pending_blocks.empty() && pending_blocks.top() == next_output_block) {
    std::string path = spillPath(next_output_block);
    {
      std::ifstream in(path, std::ios::binary);
      if (!in)
        fail("cannot open " + path);
      // Inserting an empty streambuf would set failbit on out.
      if (in.peek() != std::ifstream::traits_type::eof())
        out << in.rdbuf();
    }
    std::remove(path.c_str());
    pending_blocks.pop();
    next_output_block++;
  }
}

void processInput(int fd, StreamCursor &cursor, std::ostream &out, int width,
                  size_t &input_block_number) {
  FastReader reader;
  SDust sd;

  while (true) {
    size_t block_number;
    {
      std::lock_guard<std::mutex> lg(input_lock);
      if (!reader.LoadBlock(fd, cursor, INPUT_BLOCK_BYTES))
        break;
      block_number = input_block_number++;
    }

    reader.Parse();
    for (size_t i = 0; i < reader.size(); i++) {
      sd.seq = &reader.at(i);
      mask(sd);
      sd.reset();
    }

    // Write the block directly if it is next and the output is free;
    // otherwise spill it and let whichever thread reaches it write it.
    {
      std::unique_lock<std::mutex> lk(output_lock, std::try_to_lock);
      if (lk.owns_lock() && block_number == next_output_block) {
        flushOutput(reader, out, width);
        next_output_block++;
        writeReadySpilledBlocks(out);
        continue;
      }
    }

    std::string path = spillPath(block_number);
    {
      std::ofstream temp(path, std::ios::binary | std::ios::trunc);
      flushOutput(reader, temp, width);
      temp.close();
      if (!temp)
        fail("cannot write " + path);
    }
    std::lock_guard<std::mutex> lk(output_lock);
    pending_blocks.push(block_number);
    writeReadySpilledBlocks(out);
  }
}

int main(int argc, char **argv) {
  int ch;
  int lineWidth = 72;
  int threads = 1;
  std::string infile = "/dev/stdin";
  std::string outfile = "/dev/stdout";
  std::string buffer;
  const char *prog = "k2mask";

  struct option longopts[] = {
      {"window", required_argument, NULL, 'W'},
      {"level", required_argument, NULL, 'T'},
      {"in", required_argument, NULL, 'i'},
      {"out", required_argument, NULL, 'o'},
      {"width", required_argument, NULL, 'w'},
      {"outfmt", required_argument, NULL, 'f'},
      {"threads", required_argument, NULL, 't'},
      {"replace-masked-with", required_argument, NULL, 'r'},
      {"help", no_argument, NULL, 'h'},
      {NULL, 0, NULL, 0}};

  while ((ch = getopt_long_only(argc, argv, "W:T:hi:w:o:r:t:", longopts,
                                NULL)) != -1) {
    switch (ch) {
    case 'W':
      windowSize = atoi(optarg);
      break;
    case 'T':
      threshold = atoi(optarg);
      break;
    case 'i':
      // Input file name (default: stdin)
      infile = optarg;
      break;
    case 'w':
      // Wrap sequence every after outputting this many (72) characters
      lineWidth = atoi(optarg);
      break;
    case 'f': {
      // Output format, currently FASTA only
      std::string arg(optarg);
      std::string format("fasta");

      if (!stricasecmp(arg, format)) {
        std::cerr << prog << ":  currently only supports outputting FASTA."
                  << std::endl;
        std::exit(1);
      }
      break;
    }
    case 'o':
      // Output file name (default: stdout)
      outfile = optarg;
      break;
    case 't':
      // Number of threads. Each thread reads, masks and writes whole
      // blocks; output stays in input order.
      threads = atoi(optarg);
      break;
    case 'r': {
      // Replace masked character with a specified letter. (default tolower)
      if (strlen(optarg) != 1) {
        std::cerr << prog << ": -r expects a single character, " << optarg
                  << " given." << std::endl;
        usage(prog);
      }
      int r = optarg[0];
      processMaskedNucleotide = [=](char c) { return r; };
      break;
    }
    case 'h':
    default:
      usage(prog);
    }
  }
  argc -= optind;
  argv += optind;

  if (argc > 0) {
    std::cerr << prog << ": reads from stdin and writes to stdout" << std::endl;
    usage(prog);
  }
  auto outputFileMode = std::ios::out | std::ios::trunc;
  std::ofstream out(outfile, outputFileMode);
  if (!out)
    fail("cannot open " + outfile);

  int fd = infile == "/dev/stdin" ? fileno(stdin) : open(infile.c_str(), O_RDONLY);
  if (fd < 0)
    fail("cannot open " + infile);
  StreamCursor cursor;
  bool have_input = PrimeStream(fd, cursor);
  if (!have_input) {
    return 0;
  }

  if (threads < 1)
    threads = 1;
  if (threads > 1) {
    const char *tmpdir = getenv("TMPDIR");
    std::string tmpl = std::string(tmpdir && *tmpdir ? tmpdir : "/tmp") +
                       "/k2mask.XXXXXX";
    std::vector<char> buf(tmpl.begin(), tmpl.end());
    buf.push_back('\0');
    if (mkdtemp(buf.data()) == NULL)
      fail("cannot create a temporary directory from " + tmpl);
    spill_dir = buf.data();
  }

  // The calling thread is one of the `threads` masking threads.
  size_t input_block_number = 0;
  std::vector<std::thread> workers;
  workers.reserve(threads - 1);
  for (auto i = 0; i < threads - 1; i++) {
    workers.emplace_back(processInput, fd, std::ref(cursor), std::ref(out),
                         lineWidth, std::ref(input_block_number));
  }

  processInput(fd, cursor, out, lineWidth, input_block_number);

  for (auto &worker : workers) {
    worker.join();
  }

  if (!spill_dir.empty())
    rmdir(spill_dir.c_str());

  out.flush();
  if (!out)
    fail("cannot write " + outfile);

  return 0;
}
