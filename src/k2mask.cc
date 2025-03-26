/*
 * Symmetric Dustmasker mostly based on a similarly named implementation
 * by Heng Li as part of the minimap project.
 */

#include <algorithm>
#include <deque>
#include <functional>
#include <future>
#include <iostream>
#include <string>
#include <vector>
#include <queue>

#include <ctype.h>
#include <getopt.h>
#include <string.h>
#include <zlib.h>

#include "gzstream.h"
#include "seqreader.h"
#include "threadpool.h"

using namespace kraken2;

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
  Sequence seq;
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

void printFasta(Sequence& seq, std::ostream &out, int width) {
  out << ">";
  out.write(&seq.header[0], seq.header.size());
  if (seq.comment.size() > 0) {
    out << " ";
    out.write(&seq.comment[0], seq.comment.size());
  }
  out << '\n';
  for (size_t i = 0; i < seq.seq.size(); i += width) {
    if (i + width >= seq.seq.size())
      width = seq.seq.size() - i;
    out.write(&seq.seq[i], width);
    out << '\n';
  }
  out.flush();
}

SDust *mask(SDust *sd) {
  size_t i = 0;
  std::string &seq = sd->seq.seq;

  for (; i < seq.size(); i++) {
    if (asc2dna[(int)seq[i]] != 4) {
      int start = i;
      for (;;) {
        seq[i] = toupper(seq[i]);
        if ((i + 1) == seq.size() || asc2dna[(int)seq[i + 1]] == 4)
          break;
        i++;
      }
      char *s = &seq[0] + start;
      runSymmetricDust(*sd, s, i - start + 1, start);
      for (size_t j = 0; j < sd->ranges.size(); j++) {
        for (int i = sd->ranges[j].start; i < sd->ranges[j].finish; i++)
          s[i] = processMaskedNucleotide(s[i]);
      }
      sd->reset();
    }
  }
  return sd;
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
      // Number of threads. If _n_ threads are specified,
      // _n - 1_ threads will be used for masking and the
      // final thread used to manage I/O.
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
  std::vector<SDust *> sds(threads);
  for (size_t i = 0; i < sds.size(); i++) {
    sds[i] = new SDust();
  }
  thread_pool pool(threads - 1);
  std::queue<std::future<SDust *>> tasks;
  BatchSequenceReader reader(infile.c_str());
  for (SDust *sd = sds.back(); reader.NextSequence(sd->seq); sd = sds.back()) {
    sds.pop_back();
    if (threads > 1) {
      tasks.push(pool.submit(mask, sd));
      while (sds.empty()) {
        sd = tasks.front().get();
        printFasta(sd->seq, out, lineWidth);
        tasks.pop();
        sd->reset();
        sds.push_back(sd);
      }
    } else {
      mask(sd);
      printFasta(sd->seq, out, lineWidth);
      sds.push_back(sd);
    }
  }
  while (!tasks.empty()) {
    SDust *sd = tasks.front().get();
    printFasta(sd->seq, out, lineWidth);
    tasks.pop();
  }
  for (size_t i = 0; i < sds.size(); i++)
    delete (sds[i]);
  out.flush();
  remove("masked_sequences.txt");
}
