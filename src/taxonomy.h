/*
 * Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_TAXONOMY_H_
#define KRAKEN2_TAXONOMY_H_

#include "mmap_file.h"
#include "kraken2_headers.h"

namespace kraken2 {

// ID of a taxonomy node is implicit from its position in the node array
struct TaxonomyNode {
  uint64_t parent_id;     // Must be lower-numbered node
  uint64_t first_child;   // Must be higher-numbered node
  uint64_t child_count;   // Children of a node are in contiguous block
  uint64_t name_offset;   // Location of name in name data super-string
  uint64_t rank_offset;   // Location of rank in rank data super-string
  uint64_t external_id;   // Taxonomy ID for reporting purposes (usually NCBI)
  uint64_t godparent_id;  // Reserved for future use to enable faster traversal
};

class NCBITaxonomy {
  public:
  NCBITaxonomy(std::string nodes_filename, std::string names_filename);

  void MarkNode(uint64_t taxid);
  void ConvertToKrakenTaxonomy(const char *filename);

  private:
  std::map<uint64_t, uint64_t> parent_map_;
  std::map<uint64_t, std::string> name_map_;
  std::map<uint64_t, std::string> rank_map_;
  std::map<uint64_t, std::set<uint64_t> > child_map_;
  std::set<uint64_t> marked_nodes_;
  std::set<std::string> known_ranks_;
};

class Taxonomy {
  public:
  Taxonomy(const std::string &filename, bool memory_mapping=false);
  Taxonomy(const char *filename, bool memory_mapping=false);
  Taxonomy() : file_backed_(false), nodes_(nullptr), node_count_(0),
      name_data_(nullptr), name_data_len_(0),
      rank_data_(nullptr), rank_data_len_(0) { }
  ~Taxonomy();

  inline const TaxonomyNode *nodes() const { return nodes_; }
  inline const char *name_data() const { return name_data_; }
  inline const char *rank_data() const { return rank_data_; }
  inline const size_t node_count() const { return node_count_; }

  bool IsAAncestorOfB(uint64_t a, uint64_t b) const;
  uint64_t LowestCommonAncestor(uint64_t a, uint64_t b) const;
  void WriteToDisk(const char *filename) const;
  void MoveToMemory();

  void GenerateExternalToInternalIDMap();
  uint64_t GetInternalID(uint64_t external_id) const {
    if (external_to_internal_id_map_.count(external_id) == 0)
      return 0;
    return external_to_internal_id_map_.at(external_id);
  }

  private:
  void Init(const char *filename, bool memory_mapping);

  char const * const FILE_MAGIC = "K2TAXDAT";
  MMapFile taxonomy_data_file_;
  bool file_backed_;

  TaxonomyNode *nodes_;
  size_t node_count_;
  char *name_data_;
  size_t name_data_len_;
  char *rank_data_;
  size_t rank_data_len_;
  std::unordered_map<uint64_t, uint64_t> external_to_internal_id_map_;

  friend void NCBITaxonomy::ConvertToKrakenTaxonomy(const char *filename);
};

}

#endif  // KRAKEN_TAXONOMY_H_
