#include <map>
#include <string>

#include "kraken2_data.h"
#include "taxonomy.h"

using namespace kraken2;

extern "C" {
  Taxonomy *init_taxonomy(const char *filename) {
    Taxonomy *t = new Taxonomy(filename, false);
    t->GenerateExternalToInternalIDMap();
    return t;
  }

  taxid_t get_lca(Taxonomy *t, taxid_t a, taxid_t b) {
    auto internal_a = t->GetInternalID(a);
    auto internal_b = t->GetInternalID(b);
    auto internal_lca = t->LowestCommonAncestor(internal_a, internal_b);

    return t->nodes()[internal_lca].external_id;
  }

  taxid_t get_parent_id(Taxonomy *t, taxid_t taxid) {
    auto internal_taxid = t->GetInternalID(taxid);
    auto internal_parent = t->nodes()[internal_taxid].parent_id;

    return t->nodes()[internal_parent].external_id;
  }

  bool is_ancestor_of(Taxonomy *t, taxid_t parent, taxid_t child) {
    auto internal_parent = t->GetInternalID(parent);
    auto internal_child = t->GetInternalID(child);

    return t->IsAAncestorOfB(internal_parent, internal_child);
  }

  taxid_t get_internal_taxid(Taxonomy *t, taxid_t external_id) {
    return t->GetInternalID(external_id);
  }

  const char *get_rank(Taxonomy *t, taxid_t external_id) {
    auto internal_id = get_internal_taxid(t, external_id);
    return t->rank_data() + t->nodes()[internal_id].rank_offset;
  }

  const char *taxid_to_name(Taxonomy *t, taxid_t external_id) {
    auto internal_id = get_internal_taxid(t, external_id);
    return t->name_data() + t->nodes()[internal_id].name_offset;
  }

  taxid_t get_child_count(Taxonomy *t, taxid_t external_id) {
    auto internal_id = get_internal_taxid(t, external_id);
    return t->nodes()[internal_id].child_count;
  }

  void get_child_taxids(Taxonomy *t, taxid_t parent_taxid,
                        taxid_t *child_taxids, taxid_t num_children) {
        auto internal_id = get_internal_taxid(t, parent_taxid);
    for (auto i = 0; i < num_children; i++) {
      auto taxid = t->nodes()[internal_id].first_child + i;
      *child_taxids++ = t->nodes()[taxid].external_id;
    }
  }

  void write_to_disk(Taxonomy *t, const char *filename) {
    t->WriteToDisk(filename);
  }

  void ReadIDToTaxonMap(map<std::string, taxid_t> &id_map, std::string filename) {
    ifstream map_file(filename.c_str());
    if (! map_file.good())
      err(EX_NOINPUT, "unable to read from '%s'", filename.c_str());
    string line;

    while (getline(map_file, line)) {
      string seq_id;
      taxid_t taxid;
      istringstream iss(line);
      iss >> seq_id;
      iss >> taxid;
      if (taxid)
        id_map[seq_id] = taxid;
    }
  }

  void generate_taxonomy(const char *names, const char *nodes,
                         const char *seqid2taxid, const char *taxon_filename) {
    auto nodes_filename = std::string(nodes);
    auto names_filename = std::string(names);
    NCBITaxonomy ncbi_taxonomy(nodes_filename, names_filename);

    std::map<std::string, taxid_t> id_map;
    ReadIDToTaxonMap(id_map, seqid2taxid);
    for (auto &kv_pair : id_map) {
      if (kv_pair.second != 0) {
        ncbi_taxonomy.MarkNode(kv_pair.second);
      }
    }
    ncbi_taxonomy.ConvertToKrakenTaxonomy(taxon_filename);
  }

  void destroy_taxonomy(Taxonomy *t) {
    delete t;
  }
}
