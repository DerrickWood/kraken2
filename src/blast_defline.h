#ifndef __BLAST_DEFLINE_H__
#define __BLAST_DEFLINE_H__

#include <inttypes.h>
#include <stdio.h>

#include "kvec.h"

enum {
  seq_local,
  seq_gibbsq,
  seq_gibbmt,
  seq_giim,
  seq_genbank,
  seq_embl,
  seq_pir,
  seq_swissprot,
  seq_patent,
  seq_other,
  seq_general,
  seq_gi,
  seq_ddbj,
  seq_prf,
  seq_pdb,
  seq_tpg,
  seq_tpe,
  seq_tpd,
  seq_gpipe,
  seq_named_annot_track,
  seq_none,
};

typedef struct {
  uint8_t *buffer;
  size_t buf_len;
  size_t buf_cap;
  size_t curr_pos;
  uint32_t curr_oid;
  uint32_t *block_lens;
  uint32_t num_oids;
  FILE *asn1_file;
} asn1_t;

typedef uint32_t integer;
typedef uint8_t * visible_string;
typedef kvec_t(integer) int_vec;

typedef struct {
  visible_string name;
  visible_string acc;
  visible_string rel;
  integer ver;
} text_seqid;


typedef struct {
  integer id;
  visible_string db;
  visible_string release;
} gi_import_id;

typedef struct {
  visible_string country;
  union {
    visible_string number;
    visible_string app_number;
  } doc_number;
  int doc_number_type;
  visible_string doc_type;
} id_pat;

typedef struct {
  integer seqid;
  id_pat *cit;
} patent_seqid;

typedef struct {
  union {
    integer id;
    visible_string str;
  } id;
  int id_type;
} object_id;

typedef struct {
  integer year;
  integer month;
  integer day;
  visible_string season;
  integer hour;
  integer minute;
  integer second;
} date_std;

typedef struct {
  union {
    visible_string str;
    date_std std;
  } date_format;
  int date_type;
} date;

typedef struct {
  visible_string db;
  object_id tag;
} db_tag;

typedef struct {
  visible_string mol;
  integer chain;
  date rel;
  visible_string chain_id;
} pdb_seqid;

typedef struct {
  union {
    /* object_id local; */
    /* integer gibbsq; */
    /* integer gibbmt; */
    /* gi_import_id giim; */
    /* text_seqid genbank; */
    /* text_seqid embl; */
    /* text_seqid pir; */
    /* text_seqid swissprot; */
    /* patent_seqid patent; */
    /* text_seqid refseq; */
    /* db_tag general; */
    /* integer gi; */
    /* text_seqid ddbj; */
    /* text_seqid prf; */
    /* pdb_seqid pdb; */
    /* text_seqid tpg; */
    /* text_seqid tpe; */
    /* text_seqid tpd; */
    /* text_seqid gpipe; */
    /* text_seqid named_annot_track; */
    object_id obj_id;
    integer int_id;
    db_tag db_tag_id;
    patent_seqid pat_id;
    pdb_seqid pdb_id;
    text_seqid text_id;
    gi_import_id giim_id;
  } id;
  int seq_id_type;
} seq_id;


typedef kvec_t(seq_id) seq_id_vec;
typedef struct {
  visible_string title;
  seq_id_vec seq_ids;
  integer taxid;
  int_vec memberships;
  int_vec links;
  int_vec other_info;
} blast_defline;

typedef kvec_t(blast_defline *) blast_deflines;

asn1_t *init_asn1(const char *filename, uint32_t *block_len, uint32_t num_oids);

void free_asn1(asn1_t *asn1);

int get_asn1_visible_string(asn1_t *data, visible_string *string);

int get_asn1_integer(asn1_t *data, integer *num);

int asn1_sequence_start(asn1_t *asn1_data);

uint8_t get_explicit_tag(asn1_t *asn1_data);

int asn1_indefinite_tag_end(asn1_t *asn1_data);

int get_asn1_date(asn1_t *data, date *d);

int get_db_tag(asn1_t *data, db_tag *tag);

int get_pdb_seqid(asn1_t *data, pdb_seqid *seqid);

int get_gi_import_id(asn1_t *data, gi_import_id *seqid);

int get_text_seq_id(asn1_t *data, text_seqid *seqid);

int get_patent_seq_id(asn1_t *data, patent_seqid *seqid);

int get_object_id(asn1_t *data, object_id *id);

int get_id_pat(asn1_t *data, id_pat *id);

int get_seq_id(asn1_t *data, seq_id *seqid);

int get_ints(asn1_t *data, int_vec *vec);

int get_memberships(asn1_t *asn1_data, int_vec *memberships);

int get_links(asn1_t *asn1_data, int_vec *links);

int get_other_info(asn1_t *asn1_data, int_vec *other_info);

int get_blast_defline(asn1_t *data, blast_defline *defline);

blast_defline *init_blast_defline();

int get_blast_deflines(asn1_t *asn1_data, blast_deflines *deflines);

#endif
