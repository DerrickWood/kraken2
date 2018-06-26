/*
 * Copyright 2013-2018, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include "reports.h"

using std::vector;
using std::string;
using std::ofstream;
using std::map;

namespace kraken2 {

taxon_counts_t GetCladeCounts(Taxonomy &tax, taxon_counts_t &call_counts) {
  taxon_counts_t clade_counts;

  for (auto &kv_pair : call_counts) {
    auto taxid = kv_pair.first;
    auto count = kv_pair.second;

    while (taxid) {
      clade_counts[taxid] += count;
      taxid = tax.nodes()[taxid].parent_id;
    }
  }

  return clade_counts;
}

void PrintMpaStyleReportLine(ofstream &ofs, uint64_t clade_count, const string &taxonomy_line) {
  ofs << taxonomy_line << "\t" << clade_count << std::endl;
}

void MpaReportDFS(taxid_t taxid, ofstream &ofs, bool report_zeros, Taxonomy &taxonomy,
                  taxon_counts_t &clade_counts, vector<string> &taxonomy_names)
{
  // Clade count of 0 means all subtree nodes have clade count of 0
  if (! report_zeros && clade_counts[taxid] == 0)
    return;
  TaxonomyNode node = taxonomy.nodes()[taxid];
  string rank = taxonomy.rank_data() + node.rank_offset;

  char rank_code = '\0';
  if (rank == "superkingdom") { rank_code = 'd'; }
  else if (rank == "kingdom") { rank_code = 'k'; }
  else if (rank == "phylum")  { rank_code = 'p'; }
  else if (rank == "class")   { rank_code = 'c'; }
  else if (rank == "order")   { rank_code = 'o'; }
  else if (rank == "family")  { rank_code = 'f'; }
  else if (rank == "genus")   { rank_code = 'g'; }
  else if (rank == "species") { rank_code = 's'; }

  if (rank_code) {
    string name = "";
    name += rank_code;
    name += "__";
    name += (taxonomy.name_data() + node.name_offset);
    taxonomy_names.push_back(name);
    string taxonomy_line = "";
    for (auto &str : taxonomy_names)
      taxonomy_line += str + "|";
    taxonomy_line.pop_back();
    PrintMpaStyleReportLine(ofs, clade_counts[taxid], taxonomy_line);
  }

  auto child_count = node.child_count;
  if (child_count != 0) {
    vector<uint64_t> children(child_count);

    for (auto i = 0u; i < child_count; i++) {
      children[i] = node.first_child + i;
    }
    // Sorting child IDs by descending order of clade counts
    sort(children.begin(), children.end(),
      [&](const uint64_t &a, const uint64_t &b) {
        return clade_counts[a] > clade_counts[b];
      }
    );
    for (auto child : children) {
      MpaReportDFS(child, ofs, report_zeros, taxonomy, clade_counts,
                   taxonomy_names);
    }
  }

  if (rank_code)
    taxonomy_names.pop_back();
}

void ReportMpaStyle(string filename, bool report_zeros, Taxonomy &taxonomy, taxon_counts_t &call_counts) {
  taxon_counts_t clade_counts = GetCladeCounts(taxonomy, call_counts);
  ofstream ofs(filename);
  vector<string> taxonomy_names;
  MpaReportDFS(1, ofs, report_zeros, taxonomy, clade_counts, taxonomy_names);
}

void PrintKrakenStyleReportLine(ofstream &ofs, uint64_t total_seqs, uint64_t clade_count, uint64_t call_count,
  const string &rank_str, uint32_t taxid, const string &sci_name, int depth)
{
  char pct_buffer[10] = "";
  sprintf(pct_buffer, "%6.2f%%", 100.0 * clade_count / total_seqs);

  ofs << pct_buffer << "\t"
      << clade_count << "\t"
      << call_count << "\t"
      << rank_str << "\t"
      << taxid << "\t";
  for (auto i = 0; i < depth; i++)
    ofs << "  ";
  ofs << sci_name << std::endl;
}

// Depth-first search of taxonomy tree, reporting info at each node
void KrakenReportDFS(uint32_t taxid, ofstream &ofs, bool report_zeros,
    Taxonomy &taxonomy, taxon_counts_t &clade_counts,
    taxon_counts_t &call_counts, uint64_t total_seqs,
    char rank_code, int rank_depth, int depth)
{
  // Clade count of 0 means all subtree nodes have clade count of 0
  if (! report_zeros && clade_counts[taxid] == 0)
    return;
  TaxonomyNode node = taxonomy.nodes()[taxid];
  string rank = taxonomy.rank_data() + node.rank_offset;

  if (rank == "superkingdom") { rank_code = 'D'; rank_depth = 0; }
  else if (rank == "kingdom") { rank_code = 'K'; rank_depth = 0; }
  else if (rank == "phylum")  { rank_code = 'P'; rank_depth = 0; }
  else if (rank == "class")   { rank_code = 'C'; rank_depth = 0; }
  else if (rank == "order")   { rank_code = 'O'; rank_depth = 0; }
  else if (rank == "family")  { rank_code = 'F'; rank_depth = 0; }
  else if (rank == "genus")   { rank_code = 'G'; rank_depth = 0; }
  else if (rank == "species") { rank_code = 'S'; rank_depth = 0; }
  else { rank_depth++; }

  string rank_str(&rank_code, 0, 1);
  if (rank_depth != 0)
    rank_str += std::to_string(rank_depth);

  string name = taxonomy.name_data() + node.name_offset;

  PrintKrakenStyleReportLine(ofs, total_seqs, clade_counts[taxid],
      call_counts[taxid], rank_str, node.external_id, name, depth);

  auto child_count = node.child_count;
  if (child_count != 0) {
    vector<uint64_t> children(child_count);

    for (auto i = 0u; i < child_count; i++) {
      children[i] = node.first_child + i;
    }
    // Sorting child IDs by descending order of clade counts
    std::sort(children.begin(), children.end(),
      [&](const uint64_t &a, const uint64_t &b) {
        return clade_counts[a] > clade_counts[b];
      }
    );
    for (auto child : children) {
      KrakenReportDFS(child, ofs, report_zeros, taxonomy, clade_counts, call_counts,
          total_seqs, rank_code, rank_depth, depth + 1);
    }
  }
}

void ReportKrakenStyle(string filename, bool report_zeros, Taxonomy &taxonomy,
    taxon_counts_t &call_counts, uint64_t total_seqs, uint64_t total_unclassified)
{
  taxon_counts_t clade_counts = GetCladeCounts(taxonomy, call_counts);

  ofstream ofs(filename);

  // Special handling of the unclassified sequences
  if (total_unclassified != 0 || report_zeros)
    PrintKrakenStyleReportLine(ofs, total_seqs, total_unclassified,
                               total_unclassified, "U", 0, "unclassified", 0);
  // DFS through the taxonomy, printing nodes as encountered
  KrakenReportDFS(1, ofs, report_zeros, taxonomy, clade_counts, call_counts,
                  total_seqs, 'R', -1, 0);
}

}  // end namespace
