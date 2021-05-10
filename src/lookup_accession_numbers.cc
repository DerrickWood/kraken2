/*
 * Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include "kraken2_headers.h"
#include "mmap_file.h"
#include "utilities.h"

using namespace kraken2;
using std::string;
using std::unordered_map;
using std::vector;

int main(int argc, char **argv) {
  if (argc < 3)
    errx(EX_USAGE, "Usage: lookup_accession_numbers <lookup file> <accmaps>");

  unordered_map<string, vector<string>> target_lists;
  std::ifstream lookup_list_file(argv[1]);
  string line;
  while (getline(lookup_list_file, line)) {
    auto fields = SplitString(line, "\t", 2);
    string seqid = fields[0];
    string accnum = fields[1];
    target_lists[accnum].push_back(seqid);
  }
  lookup_list_file.close();

  auto initial_target_count = target_lists.size();

  MMapFile accmap_file;
  if (isatty(fileno(stderr)))
    std::cerr << "\rFound 0/" << initial_target_count << " targets...";
  uint64_t accessions_searched = 0;
  for (int i = 2; i < argc; i++) {
    if (target_lists.empty())  // Stop processing files if we've found all we need
      break;
    accmap_file.OpenFile(argv[i]);
    char *ptr = accmap_file.fptr();
    // Skip header line
    char *lf_ptr = (char *) memchr(ptr, '\n', accmap_file.filesize());
    if (lf_ptr != nullptr)
      ptr = lf_ptr + 1;

    while ((size_t)(ptr - accmap_file.fptr()) < accmap_file.filesize()) {
      lf_ptr = (char *) memchr(ptr, '\n', accmap_file.filesize() - (ptr - accmap_file.fptr()));
      if (lf_ptr == nullptr) {
        warnx("expected EOL not found at EOF in %s", argv[i]);
        break;
      }
      char *tab_ptr = (char *) memchr(ptr, '\t', lf_ptr - ptr);
      if (tab_ptr == nullptr) {
        warnx("expected TAB not found in %s", argv[i]);
        break;
      }
      string accnum(ptr, tab_ptr - ptr);
      accessions_searched++;
      if (target_lists.count(accnum)) {
        for (int j = 0; j < 2; j++) {
          ptr = tab_ptr + 1;
          tab_ptr = (char *) memchr(ptr, '\t', lf_ptr - ptr);
          if (tab_ptr == nullptr) {
            warnx("expected TAB not found in %s", argv[i]);
            break;
          }
        }
        if (tab_ptr == nullptr)
          break;
        string taxid(ptr, tab_ptr - ptr);
        for (auto &seqid : target_lists[accnum])
          std::cout << seqid << "\t" << taxid << std::endl;
        target_lists.erase(accnum);
        if (isatty(fileno(stderr)))
          std::cerr << "\rFound " << (initial_target_count - target_lists.size())
              << "/" << initial_target_count << " targets, searched through "
              << accessions_searched << " accession IDs...";
        if (target_lists.empty())  // Stop processing file if we've found all we need
          break;
      }
      if (accessions_searched % 10000000 == 0 && isatty(fileno(stderr)))
        std::cerr << "\rFound " << (initial_target_count - target_lists.size())
            << "/" << initial_target_count << " targets, searched through "
            << accessions_searched << " accession IDs...";
      ptr = lf_ptr + 1;
    }
    accmap_file.CloseFile();
  }
  if (isatty(fileno(stderr)))
    std::cerr << "\r";
  std::cerr << "Found " << (initial_target_count - target_lists.size())
      << "/" << initial_target_count << " targets, searched through "
      << accessions_searched << " accession IDs, search complete." << std::endl;

  if (! target_lists.empty()) {
    std::cerr << "lookup_accession_numbers: " << target_lists.size() << "/"
         << initial_target_count << " accession numbers remain unmapped, see "
         << "unmapped.txt in DB directory" << std::endl;
    std::ofstream ofs("unmapped.txt");
    for (auto &kv_pair : target_lists)
      ofs << kv_pair.first << std::endl;
  }

  return 0;
}
