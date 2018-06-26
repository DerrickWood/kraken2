/*
 * Copyright 2013-2018, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include "aa_translate.h"

using std::string;
using std::vector;

namespace kraken2 {

/*
This *was* the map, which was ridiculously slow.  Note that ordering is
AGCT, not ACGT.
static map<string, char> old_translation_map = {
  { "AAA", 'K' }, { "AAG", 'K' }, { "AAC", 'N' }, { "AAT", 'N' },
  { "AGA", 'R' }, { "AGG", 'R' }, { "AGC", 'S' }, { "AGT", 'S' },
  { "ACA", 'T' }, { "ACG", 'T' }, { "ACC", 'T' }, { "ACT", 'T' },
  { "ATA", 'I' }, { "ATG", 'M' }, { "ATC", 'I' }, { "ATT", 'I' },
  { "GAA", 'E' }, { "GAG", 'E' }, { "GAC", 'D' }, { "GAT", 'D' },
  { "GGA", 'G' }, { "GGG", 'G' }, { "GGC", 'G' }, { "GGT", 'G' },
  { "GCA", 'A' }, { "GCG", 'A' }, { "GCC", 'A' }, { "GCT", 'A' },
  { "GTA", 'V' }, { "GTG", 'V' }, { "GTC", 'V' }, { "GTT", 'V' },
  { "CAA", 'Q' }, { "CAG", 'Q' }, { "CAC", 'H' }, { "CAT", 'H' },
  { "CGA", 'R' }, { "CGG", 'R' }, { "CGC", 'R' }, { "CGT", 'R' },
  { "CCA", 'P' }, { "CCG", 'P' }, { "CCC", 'P' }, { "CCT", 'P' },
  { "CTA", 'L' }, { "CTG", 'L' }, { "CTC", 'L' }, { "CTT", 'L' },
  { "TAA", '*' }, { "TAG", '*' }, { "TAC", 'Y' }, { "TAT", 'Y' },
  { "TGA", '*' }, { "TGG", 'W' }, { "TGC", 'C' }, { "TGT", 'C' },
  { "TCA", 'S' }, { "TCG", 'S' }, { "TCC", 'S' }, { "TCT", 'S' },
  { "TTA", 'L' }, { "TTG", 'L' }, { "TTC", 'F' }, { "TTT", 'F' }
};
*/

static char translation_map[] = "KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLL**YY*WCCSSSSLLFF";

void TranslateToAllFrames(string &dna_seq, vector<string> &aa_seqs) {
  auto max_size = (dna_seq.size() / 3) + 1;
  for (auto i = 0; i < 6; i++)
    aa_seqs[i].assign(max_size, ' ');
  if (dna_seq.size() < 3)
    return;

  uint8_t fwd_codon = 0, rev_codon = 0;
  int ambig_nt_countdown = 0;  // if positive, bases to go until N leaves codon
  size_t frame_len[6] = {0};
  for (auto i = 0u; i < dna_seq.size(); i++) {
    auto frame = i % 3;
    fwd_codon <<= 2;
    fwd_codon &= 0x3f;
    rev_codon >>= 2;
    if (ambig_nt_countdown)
      ambig_nt_countdown--;
    // Map is based on AGCT coding, not ACGT
    switch (dna_seq[i]) {
      case 'A' :                    rev_codon |= 0x30; break;
      case 'G' : fwd_codon |= 0x01; rev_codon |= 0x20; break;
      case 'C' : fwd_codon |= 0x02; rev_codon |= 0x10; break;
      case 'T' : fwd_codon |= 0x03; break;
      default:
        ambig_nt_countdown = 3;
    }

    if (i >= 2) {  // we've got a full codon
      char ch;

      // Translate and append to forward frame
      ch = ambig_nt_countdown == 0 ? translation_map[fwd_codon] : 'X';
      aa_seqs[frame][frame_len[frame]++] = ch;

      // Translate and prepend to reverse frame
      ch = ambig_nt_countdown == 0 ? translation_map[rev_codon] : 'X';
      aa_seqs[frame + 3][max_size - 1 - frame_len[frame + 3]++] = ch;
    }
  }

  for (auto i = 0; i < 3; i++)
    if (aa_seqs[i].size() != frame_len[i])
      aa_seqs[i].resize(frame_len[i]);
  for (auto i = 3; i < 6; i++)
    if (aa_seqs[i].size() != frame_len[i])
      aa_seqs[i].erase(0, aa_seqs[i].size() - frame_len[i]);
}

}
