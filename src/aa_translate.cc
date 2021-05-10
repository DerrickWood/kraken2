/*
 * Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
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
static uint8_t fwd_lookup_table[UINT8_MAX + 1] = {0};
static uint8_t rev_lookup_table[UINT8_MAX + 1] = {0};

void TranslateToAllFrames(string &dna_seq, vector<string> &aa_seqs) {
  auto max_size = (dna_seq.size() / 3) + 1;
  for (auto i = 0; i < 6; i++)
    aa_seqs[i].assign(max_size, ' ');
  if (dna_seq.size() < 3)
    return;

  if (fwd_lookup_table[0] == 0) {
    for (size_t i = 0; i <= UINT8_MAX; i++)
      fwd_lookup_table[i] = rev_lookup_table[i] = UINT8_MAX;
    // Map is based on AGCT coding, not ACGT
    fwd_lookup_table[(int) 'A'] = fwd_lookup_table[(int) 'a'] = 0x00;
    fwd_lookup_table[(int) 'G'] = fwd_lookup_table[(int) 'g'] = 0x01;
    fwd_lookup_table[(int) 'C'] = fwd_lookup_table[(int) 'c'] = 0x02;
    fwd_lookup_table[(int) 'T'] = fwd_lookup_table[(int) 't'] = 0x03;
    rev_lookup_table[(int) 'A'] = rev_lookup_table[(int) 'a'] = 0x30;
    rev_lookup_table[(int) 'G'] = rev_lookup_table[(int) 'g'] = 0x20;
    rev_lookup_table[(int) 'C'] = rev_lookup_table[(int) 'c'] = 0x10;
    rev_lookup_table[(int) 'T'] = rev_lookup_table[(int) 't'] = 0x00;
  }

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
    auto fwd_lookup_code = fwd_lookup_table[(int) dna_seq[i]];
    auto rev_lookup_code = rev_lookup_table[(int) dna_seq[i]];
    if (fwd_lookup_code == UINT8_MAX)
      ambig_nt_countdown = 3;
    else {
      fwd_codon |= fwd_lookup_code;
      rev_codon |= rev_lookup_code;
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
