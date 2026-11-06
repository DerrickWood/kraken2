// Exercise the real production function without adding timing to src/.
// Unused classifier entry points are discarded by the benchmark link command.
#define main kraken2_classifier_main
#include "../../src/classify.cc"
#undef main

#include <chrono>
#include <iomanip>
#include <random>
#include <stdexcept>

#include "resolve_tree_reference.h"

namespace {
using Clock = std::chrono::steady_clock;
static_assert(Clock::is_steady, "benchmark requires a monotonic clock");
using Resolver = taxid_t (*)(taxon_counts_t &, Taxonomy &, size_t, const Options &);
volatile uint64_t checksum_sink = 0;
size_t checked = 0;

struct Input {
  taxon_counts_t hits;
  size_t total = 0;
};

// Use the existing on-disk taxonomy layout; no private access or mock LCA.
void WriteTree(const std::string &path, const std::vector<uint64_t> &parents) {
  std::vector<TaxonomyNode> nodes(parents.size());
  for (size_t i = 1; i < nodes.size(); ++i) {
    nodes[i].parent_id = parents[i];
    nodes[i].external_id = i;
  }
  std::ofstream out(path, std::ios::binary);
  const size_t count = nodes.size(), zero = 0;
  out.write("K2TAXDAT", 8);
  out.write(reinterpret_cast<const char *>(&count), sizeof(count));
  out.write(reinterpret_cast<const char *>(&zero), sizeof(zero));
  out.write(reinterpret_cast<const char *>(&zero), sizeof(zero));
  out.write(reinterpret_cast<const char *>(nodes.data()), nodes.size() * sizeof(nodes[0]));
  if (!out)
    throw std::runtime_error("cannot write taxonomy fixture");
}

void Check(const Input &input, Taxonomy &tax, double confidence) {
  auto before = input.hits, after = input.hits;
  Options opts{};
  opts.confidence_threshold = confidence;
  const auto expected = ResolveTreeReference(before, tax, input.total, opts);
  const auto actual = ResolveTree(after, tax, input.total, opts);
  if (expected != actual || before != after) {
    std::ostringstream message;
    message << "ResolveTree mismatch: D=" << input.hits.size()
            << " confidence=" << confidence << " expected=" << expected
            << " actual=" << actual;
    throw std::runtime_error(message.str());
  }
  ++checked;
}

Input RandomInput(size_t distinct, size_t nodes, std::mt19937_64 &rng) {
  Input input;
  while (input.hits.size() < distinct)
    input.hits.emplace(1 + rng() % (nodes - 1), 1 + rng() % 20);
  for (const auto &hit : input.hits)
    input.total += hit.second;
  input.total += rng() % (input.total + 1);  // unmatched positions in denominator
  return input;
}

void CheckTree(Taxonomy &tax, std::mt19937_64 &rng) {
  const double thresholds[] = {0, 0.01, 0.5, 0.9, 1};
  for (double threshold : thresholds) {
    Check({}, tax, threshold);
    Check({{{1, 1}}, 1}, tax, threshold);
    Check({{{0, 0}, {2, 0}, {3, 0}}, 10}, tax, threshold);
    Check({{{2, 10}, {3, 10}}, 20}, tax, threshold);
    Check({{{1, 3}, {2, 8}, {3, 8}}, 30}, tax, threshold);
    // Keep the original uint32_t score wraparound semantics as well.
    Check({{{1, UINT64_C(0xffffffff)}, {2, 2}, {3, UINT64_C(0x100000005)}},
           1000}, tax, threshold);
  }
  for (size_t i = 0; i < 1200; ++i) {
    auto input = RandomInput(rng() % 33, tax.node_count(), rng);
    if (i % 7 == 0)
      input.hits[0] = 0;
    if (i % 11 == 0)
      for (auto &hit : input.hits)
        hit.second = rng();  // modulo-2^32 scoring, including zero sentinel
    const double threshold = thresholds[i % 5];
    Check(input, tax, threshold);
    // Reinsert in another order to exercise ties and hash-map iteration order.
    std::vector<std::pair<taxid_t, uint64_t>> entries(input.hits.begin(), input.hits.end());
    std::shuffle(entries.begin(), entries.end(), rng);
    input.hits.clear();
    input.hits.rehash(97);
    for (const auto &entry : entries)
      input.hits.insert(entry);
    Check(input, tax, threshold);
  }
}

// A batch timer keeps clock overhead outside individual ResolveTree calls.
// Input construction, map copying, correctness checks and I/O are not timed.
double Measure(Resolver fn, std::vector<Input> &inputs, Taxonomy &tax,
               const Options &opts, size_t repeats) {
  Resolver volatile invoke = fn;
  uint64_t checksum = 0;
  const auto start = Clock::now();
  for (size_t r = 0; r < repeats; ++r)
    for (auto &input : inputs)
      checksum += invoke(input.hits, tax, input.total, opts);
  const auto stop = Clock::now();
  checksum_sink = checksum;
  return std::chrono::duration<double, std::nano>(stop - start).count();
}

size_t Calibrate(Resolver fn, std::vector<Input> &inputs, Taxonomy &tax,
                 const Options &opts) {
  size_t repeats = 1;
  while (Measure(fn, inputs, tax, opts, repeats) < 5000000 && repeats < 65536)
    repeats *= 2;
  return repeats;
}

double Median(std::vector<double> values) {
  std::sort(values.begin(), values.end());
  return values[values.size() / 2];
}

void Benchmark(const std::string &shape, size_t distinct, double threshold,
               Taxonomy &tax, std::mt19937_64 &rng, size_t rounds) {
  std::vector<Input> before, after;
  for (size_t i = 0; i < 16; ++i) {
    before.push_back(RandomInput(distinct, tax.node_count(), rng));
    Check(before.back(), tax, threshold);
  }
  after = before;
  Options opts{};
  opts.confidence_threshold = threshold;
  // Warm both functions, including the legacy hit_counts[max_taxon] insertion.
  Measure(ResolveTreeReference, before, tax, opts, 1);
  Measure(ResolveTree, after, tax, opts, 1);
  const auto old_repeats = Calibrate(ResolveTreeReference, before, tax, opts);
  const auto new_repeats = Calibrate(ResolveTree, after, tax, opts);
  std::vector<double> old_times, new_times;
  for (size_t r = 0; r < rounds; ++r) {
    auto old_sample = [&] {
      old_times.push_back(Measure(ResolveTreeReference, before, tax, opts, old_repeats)
                          / (old_repeats * before.size()));
    };
    auto new_sample = [&] {
      new_times.push_back(Measure(ResolveTree, after, tax, opts, new_repeats)
                          / (new_repeats * after.size()));
    };
    if (r % 2) { new_sample(); old_sample(); }
    else { old_sample(); new_sample(); }
  }
  const auto old_ns = Median(old_times), new_ns = Median(new_times);
  std::cout << shape << ',' << distinct << ',' << threshold << ','
            << old_ns << ',' << new_ns << ',' << old_ns / new_ns << '\n';
}
}  // namespace

int main(int argc, char **argv) {
  try {
    if (argc != 4)
      throw std::runtime_error("usage: benchmark TEMP_DIR ROUNDS CHECK_ONLY");
    const size_t rounds = std::stoul(argv[2]);
    const bool check_only = std::string(argv[3]) == "1";
    if (rounds == 0)
      throw std::runtime_error("rounds must be positive");
    std::mt19937_64 rng(20261022);
    std::cout << "shape,distinct_taxa,confidence,reference_ns,current_ns,speedup\n"
              << std::fixed << std::setprecision(3);
    for (const std::string shape : {"balanced", "chain", "star", "irregular"}) {
      const size_t count = shape == "chain" ? 257 : 4097;
      std::vector<uint64_t> parents(count);
      for (size_t i = 2; i < count; ++i) {
        if (shape == "chain") parents[i] = i - 1;
        else if (shape == "star") parents[i] = 1;
        else if (shape == "balanced") parents[i] = (i - 2) / 4 + 1;
        else parents[i] = 1 + rng() % (i - 1);
      }
      const auto path = std::string(argv[1]) + "/" + shape + ".k2d";
      WriteTree(path, parents);
      Taxonomy tax(path, false);
      CheckTree(tax, rng);
      if (!check_only)
        for (size_t distinct : {1, 2, 4, 8, 16, 64, 256})
          for (double threshold : {0.0, 0.5, 0.9})
            Benchmark(shape, distinct, threshold, tax, rng, rounds);
    }
    // Exhaust all 0/1/2 vote assignments on a small tree, with missing keys.
    const auto path = std::string(argv[1]) + "/exhaustive.k2d";
    WriteTree(path, {0, 0, 1, 1, 2, 2, 3, 3});
    Taxonomy tax(path, false);
    for (size_t assignment = 0; assignment < 2187; ++assignment) {
      Input input;
      size_t digits = assignment;
      for (taxid_t id = 1; id <= 7; ++id, digits /= 3) {
        const auto votes = digits % 3;
        if (votes) input.hits[id] = votes;
        input.total += votes;
      }
      for (double threshold : {0.0, 0.01, 0.5, 0.9, 1.0}) {
        Check(input, tax, threshold);
        auto unmatched = input;
        ++unmatched.total;
        Check(unmatched, tax, threshold);
      }
    }
    std::cout << "# correctness: " << checked << " cases passed (taxid and map contents)\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
