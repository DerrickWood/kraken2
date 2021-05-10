/*
 * Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_REPORTS_H_
#define KRAKEN2_REPORTS_H_

#include <unordered_map>
#include <iomanip>
#include "kraken2_headers.h"
#include "taxonomy.h"
#include "kraken2_data.h"
#include "readcounts.h"

namespace kraken2 {

// Still TODO: Create an MPA-style reporter that can take a std::vector of
//   call_counts and a std::vector of sample IDs
void ReportMpaStyle(std::string filename, bool report_zeros, Taxonomy &tax, taxon_counters_t &call_counters);
taxon_counts_t GetCladeCounts(Taxonomy &tax, taxon_counts_t &call_counts);
taxon_counters_t GetCladeCounters(Taxonomy &tax, taxon_counters_t &call_counters);
void PrintMpaStyleReportLine(std::ofstream &ofs, uint64_t clade_count, const std::string &taxonomy_line);
void MpaReportDFS(taxid_t taxid, std::ofstream &ofs, bool report_zeros, Taxonomy &taxonomy,
    taxon_counts_t &clade_counts, std::vector<std::string> &taxonomy_names);
void PrintKrakenStyleReportLine(std::ofstream &ofs, bool report_kmer_data,
    uint64_t total_seqs, READCOUNTER clade_counter, READCOUNTER taxon_counter,
    const std::string &rank_str, uint32_t taxid, const std::string &sci_name, int depth);
void KrakenReportDFS(uint32_t taxid, std::ofstream &ofs, bool report_zeros,
    bool report_kmer_data, Taxonomy &taxonomy, taxon_counters_t &clade_counters,
    taxon_counters_t &call_counters, uint64_t total_seqs, char rank_code, int rank_depth, int depth);
void ReportKrakenStyle(std::string filename, bool report_zeros, bool report_kmer_data,
    Taxonomy &taxonomy, taxon_counters_t &call_counters, uint64_t total_seqs, uint64_t total_unclassified);

}

#endif
