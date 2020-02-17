/*
 * Copyright 2013-2020, Derrick Wood <dwood@cs.jhu.edu>
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
void ReportMpaStyle(std::string filename, bool report_zeros, Taxonomy &tax, taxon_counts_t &call_counts);
//void ReportKrakenStyle(std::string filename, Taxonomy &tax, taxon_counts_t &call_counts, uint64_t total_unclassified);
taxon_counts_t GetCladeCounts(Taxonomy &tax, taxon_counts_t &call_counts);
void PrintMpaStyleReportLine(std::ofstream &ofs, uint64_t clade_count, const std::string &taxonomy_line);
void MpaReportDFS(taxid_t taxid, std::ofstream &ofs, bool report_zeros, Taxonomy &taxonomy,
    taxon_counts_t &clade_counts, std::vector<std::string> &taxonomy_names);
void PrintKrakenStyleReportLine(std::ofstream &ofs, uint64_t total_seqs, uint64_t clade_count, uint64_t call_count,
    const std::string &rank_str, uint32_t taxid, const std::string &sci_name, int depth);
void KrakenReportDFS(uint32_t taxid, std::ofstream &ofs, bool report_zeros,
    Taxonomy &taxonomy, taxon_counts_t &clade_counts,
    uint64_t total_seqs, char rank_code, int rank_depth, int depth);
void ReportKrakenStyle(std::string filename, bool report_zeros, Taxonomy &taxonomy,
    taxon_counts_t &call_counts, 
    uint64_t total_seqs, uint64_t total_unclassified);

    template<typename READCOUNTS>
    class TaxReport {
    private:
        const Taxonomy * taxonomy;
        const unordered_map<taxid_t, READCOUNTS>& _taxCounts; // set in constructor, from classification
        const bool _show_zeros;
        unordered_map<taxid_t, vector<const READCOUNTS*> > _children;
        unordered_map<taxid_t, READCOUNTS> _cladeCounts; // consider accessing by TaxEntry*
    public:

        TaxReport(const Taxonomy & taxonomy,
                  const unordered_map<taxid_t, READCOUNTS>& readCounts,
                  bool show_zeros)
                : taxonomy(&taxonomy),  _taxCounts(readCounts), _show_zeros(show_zeros) {
            cerr << "Setting values in the taxonomy tree ...";
            //unordered_map<const TaxonomyEntry<TAXID>*, unordered_set<const TaxonomyEntry<TAXID>*> > _children1;
            for (auto it = _taxCounts.begin(); it != _taxCounts.end(); ++it) {
                TaxonomyNode node;
                taxid_t tax = it->first;
                while (tax != 0) {
                    node = taxonomy.nodes()[tax];
                    _children[tax].push_back(&(it->second));
                    tax = node.parent_id;
                }
            }

#pragma omp parallel for schedule(dynamic, 50)
            for (size_t i = 0; i < _children.size(); ++i) {

                auto cit = _children.begin();
                advance(cit, i);
                READCOUNTS rc = *(cit->second.front());
                for (size_t j = 1; j < cit->second.size(); ++j){
                    auto counter = cit->second[j];
                    rc += *counter;
                }

#pragma omp critical(update_clade_counts)
                {
                    _cladeCounts.insert( std::make_pair( cit->first, std::move(rc) ) );
                }
            }

            cerr << " done" << endl;
        }

        void ReportKrakenStyle(std::string filename, uint64_t total_sequences, uint64_t total_unclassified){
            ofstream ofs(filename);

//            // Special handling of the unclassified sequences
            if (total_unclassified != 0 || _show_zeros) {
//                printLine(ofstream& ofs, const taxid_t taxid,
//                               const READCOUNTS & rc,
//                               unsigned depth, const string &rank_str, uint32_t external_taxid,
//                               const string &sci_name, uint64_t total_sequences)
                READCOUNTS rc(total_unclassified, 0);
                printLine(ofs, 0, rc,  0, "U", 0, "unclassified", total_sequences);
            }

            printReport(ofs, 1, 0,  'R', -1, total_sequences);
        }

        void printReport(ofstream & ofs, const taxid_t taxid, unsigned depth, char rank_code, int rank_depth, uint64_t total_sequences) {
            const auto taxit_ptr = _cladeCounts.find(taxid);
            if (taxit_ptr == _cladeCounts.end())
                return;
            const auto & cladecounts = taxit_ptr->second;
            if (!_show_zeros && cladecounts.readCount() == 0)
                return;

            TaxonomyNode node = taxonomy->nodes()[taxid];
            string rank = taxonomy->rank_data() + node.rank_offset;

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

            string name = taxonomy->name_data() + node.name_offset;

            printLine(ofs, taxid, cladecounts, depth, rank_str, node.external_id, name, total_sequences);


            // Sort children
            vector<size_t> pos;
            unordered_map<size_t, READCOUNTS *> rc;
            auto child_count = node.child_count;
            if (child_count != 0) {
                vector<uint64_t> children(child_count);

                for (auto i = 0u; i < child_count; i++) {
                    children[i] = node.first_child + i;
                }

                for (size_t i = 0; i < children.size(); ++i) {
                    auto it = _cladeCounts.find(children[i]);
                    if (it != _cladeCounts.end()) {
                        pos.push_back(i);
                        rc[i] = &(it->second);
                    }
                }
                std::sort(pos.begin(), pos.end(), [&](size_t a, size_t b) { return *(rc.at(b)) < *(rc.at(a)); });

                for (size_t i = 0; i < rc.size(); ++i) {
                    auto child_taxid = children[pos[i]];
                    printReport(ofs, child_taxid, depth + 1, rank_code, rank_depth, total_sequences);
                }
            }
        }

        void printLine(ofstream& ofs, const taxid_t taxid,
                       const READCOUNTS & rc,
                       unsigned depth, const string &rank_str, uint32_t external_taxid,
                       const string &sci_name, uint64_t total_sequences)
        {
            const auto r_it = _taxCounts.find(taxid);
            const bool has_tax_data = r_it != _taxCounts.end();

            ofs << setprecision(4) << 100.0*rc.readCount()/total_sequences << "\t"
                << rc.readCount() << "\t"
                << (has_tax_data? r_it->second.readCount() : 0) << "\t"
                << rc.kmerCount() << "\t"
                << rc.uniqueKmerCount() << "\t"
                << rank_str << "\t"
                << external_taxid << "\t" ;
            for (auto i = 0; i < depth; i++)
                ofs << "  ";
            ofs << sci_name << std::endl;
        }
    };


}

#endif
