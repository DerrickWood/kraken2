#ifndef __BLAST_TO_FASTA_H__
#define __BLAST_TO_FASTA_H__

#include <inttypes.h>
#include <stdio.h>

#include "blast_defline.h"

typedef union {
  struct {
    uint8_t tag : 5;
    uint8_t construction : 1;
    uint8_t tag_class : 2;

  };
  uint8_t raw_tag;
} asn1_tag;

typedef struct {
  uint64_t offset : 32;
  uint64_t unused : 16;
  uint64_t length : 12;
  uint64_t value : 4;
} amb64_t;

typedef struct {
  uint32_t offset : 24;
  uint32_t length : 4;
  uint32_t value : 4;
} amb32_t;

typedef struct {
  union {
    amb32_t amb32;
    amb64_t amb64;
  };
} amb_t;

typedef struct {
  uint8_t *string;
  uint32_t cap;
  uint32_t len;
} string_t;

typedef struct {
  FILE *idx_file;
  uint32_t fmt_version;
  uint32_t db_seq_type;
  uint32_t volume;
  string_t title;
  string_t lmdb_file;
  string_t date;
  uint32_t num_oids;
  uint64_t vol_len;
  uint32_t max_seq_len;
  uint32_t *hdr_arr;
  uint32_t *seq_arr;
  uint32_t *amb_arr;
} idx_t;

typedef struct {
  FILE *hdr_file;
  string_t fasta_hdr;
  blast_deflines deflines;
  asn1_t *asn1;
} hdr_t;

typedef struct {
  FILE *seq_file;
  string_t buffer;
  string_t seq;
  amb_t *amb_data;
  uint32_t curr_pos;
  uint32_t amb_data_cap;
} seq_t;

#endif
