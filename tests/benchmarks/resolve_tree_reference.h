// Frozen pre-optimization ResolveTree from upstream commit 8e35df55896c672b340a99b05f55cc809db2c5ac.
// Copyright 2013-2023, Derrick Wood; distributed under the repository license.
// Keep this independent reference unchanged when optimizing the production function.

taxid_t ResolveTreeReference(taxon_counts_t &hit_counts,
    Taxonomy &taxonomy, size_t total_minimizers, const Options &opts)
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

