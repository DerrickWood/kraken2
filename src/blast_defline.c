#include <assert.h>
#include <err.h>
#include <strings.h>

#include "blast_defline.h"
#include "blast_utils.h"

asn1_t *init_asn1(const char *filename, uint32_t *block_lens, uint32_t num_oids) {
  asn1_t *asn1 = alloc_memory(NULL, sizeof(asn1_t), 0, 1);
  asn1->asn1_file = open_file(filename, "r");
  asn1->buffer = NULL;
  asn1->buf_len = 0;
  asn1->buf_cap = 0;
  asn1->curr_pos = 0;
  asn1->curr_oid = 0;
  asn1->block_lens = block_lens;
  asn1->num_oids = num_oids;

  return asn1;
}

void free_asn1(asn1_t *asn1) {
  if (asn1 != NULL) {
    fclose(asn1->asn1_file);
    if (asn1->buffer != NULL) {
      free(asn1->buffer);
    }
    free(asn1);
  }
}

size_t load_more_data(asn1_t *asn1_data) {
  if (asn1_data->curr_oid == asn1_data->num_oids) {
    return 0;
  }

  uint32_t block_len = asn1_data->block_lens[asn1_data->curr_oid++];
  if ((asn1_data->buf_cap - asn1_data->buf_len) < block_len) {
    size_t new_cap = next_power_of_2(asn1_data->buf_len + block_len);
    asn1_data->buffer = alloc_memory(asn1_data->buffer, 1, asn1_data->buf_cap, new_cap);
    asn1_data->buf_cap = new_cap;
  }
  read_into_buffer(asn1_data->asn1_file, asn1_data->buffer + asn1_data->buf_len,
                   1, block_len);
  asn1_data->buf_len += block_len;

  return block_len;
}
uint8_t asn1_get_byte(asn1_t *asn1_data) {
  // add some error checking here if the buffer is exhausted
  if (asn1_data->curr_pos == asn1_data->buf_len) {
    load_more_data(asn1_data);
  }

  return asn1_data->buffer[asn1_data->curr_pos++];
}

void asn1_backtrack(asn1_t *asn1_data) {
  asn1_data->curr_pos--;
}

uint8_t asn1_peek_byte(asn1_t *asn1_data) {
  // add some error checking here if the buffer is exhausted
  if (asn1_data->curr_pos == asn1_data->buf_len) {
    load_more_data(asn1_data);
  }

  return asn1_data->buffer[asn1_data->curr_pos];
}

int asn1_sequence_start(asn1_t *asn1_data) {
  if (asn1_peek_byte(asn1_data) != 0x30) {
    return 0;
  }

  uint8_t tag = asn1_get_byte(asn1_data);
  uint8_t len = asn1_get_byte(asn1_data);

  return len == 0x80;
}

int asn1_indefinite_tag_end(asn1_t *asn1_data) {
  uint8_t null_octet1 = asn1_get_byte(asn1_data);
  uint8_t null_octet2 = asn1_get_byte(asn1_data);

  return null_octet1 == 0 && null_octet2 == 0;
}

int asn1_is_end_of_sequence(asn1_t *asn1_data) {
  int end = asn1_indefinite_tag_end(asn1_data);
  asn1_backtrack(asn1_data);
  asn1_backtrack(asn1_data);

  return end;
}

int asn1_get_visible_string(asn1_t *asn1_data, visible_string *string) {
  uint32_t str_len = 0;
  uint8_t tag = asn1_get_byte(asn1_data);
  assert(tag == 0x1A);
  uint32_t octets = asn1_get_byte(asn1_data) & 0xff;

  /* The string length spans multiple octets */
  if ((octets & 0x80) == 0x80) {
    uint32_t num_octets = octets & (0x80 - 1);
    for (int j = 0; j < num_octets; j++) {
      octets = asn1_get_byte(asn1_data);
      str_len = (str_len << (j * 8)) | (octets & 0xff);
    }
  } else {
    str_len = octets;
  }
  if ((asn1_data->curr_pos + str_len) > asn1_data->buf_len) {
    load_more_data(asn1_data);
  }

  *string = asn1_data->buffer + asn1_data->curr_pos;
  asn1_data->curr_pos += str_len;

  return str_len;
}

int asn1_get_integer(asn1_t *asn1_data, integer *num) {
  uint32_t value = 0;
  uint8_t tag = asn1_get_byte(asn1_data);
  assert(tag = 0x02);

  uint32_t int_len = asn1_get_byte(asn1_data);
  value = asn1_get_byte(asn1_data);
  for (uint32_t j = 1; j < int_len; j++) {
    value = (value << 8) + asn1_get_byte(asn1_data);
  }

  *num = value;

  return 1;
}

uint8_t get_explicit_tag(asn1_t *asn1_data) {
  uint8_t field_no = asn1_get_byte(asn1_data);
  (void)asn1_get_byte(asn1_data);

  return field_no;
}

void asn1_get_optional_visible_string_field(asn1_t *asn1_data, uint8_t field_tag,
                                            visible_string *destination) {
  uint8_t tag = asn1_peek_byte(asn1_data);
  if (tag == field_tag) {
    (void)get_explicit_tag(asn1_data);
    asn1_get_visible_string(asn1_data, destination);
    asn1_indefinite_tag_end(asn1_data);
  }
}

void asn1_get_optional_integer_field(asn1_t *asn1_data, uint8_t field_tag,
                                            integer *destination) {
  uint8_t tag = asn1_peek_byte(asn1_data);
  if (tag == field_tag) {
    (void)get_explicit_tag(asn1_data);
    asn1_get_integer(asn1_data, destination);
    asn1_indefinite_tag_end(asn1_data);
  }
}

void asn1_get_mandatory_integer_field(asn1_t *asn1_data, integer *destination) {
    (void)get_explicit_tag(asn1_data);
    asn1_get_integer(asn1_data, destination);
    asn1_indefinite_tag_end(asn1_data);
}

void asn1_get_mandatory_visible_string_field(asn1_t *asn1_data, visible_string *destination) {
    (void)get_explicit_tag(asn1_data);
    asn1_get_visible_string(asn1_data, destination);
    asn1_indefinite_tag_end(asn1_data);
}

int get_date(asn1_t *asn1_data, date *d) {
  uint8_t field_no = asn1_peek_byte(asn1_data);
  if (field_no == 0xA0) {
    asn1_get_mandatory_visible_string_field(asn1_data, &d->date_format.str);
    d->date_type = 1;
  } else {
    (void)get_explicit_tag(asn1_data);
    (void)asn1_sequence_start(asn1_data);
    asn1_get_mandatory_integer_field(asn1_data, &d->date_format.std.year);
    asn1_get_optional_integer_field(asn1_data, 0xA1, &d->date_format.std.month);
    asn1_get_optional_integer_field(asn1_data, 0xA2, &d->date_format.std.day);
    asn1_get_optional_visible_string_field(asn1_data, 0xA3,
                                           &d->date_format.std.season);
    asn1_get_optional_integer_field(asn1_data, 0xA4, &d->date_format.std.hour);
    asn1_get_optional_integer_field(asn1_data, 0xA5,
                                    &d->date_format.std.minute);
    asn1_get_optional_integer_field(asn1_data, 0xA6,
                                    &d->date_format.std.second);
    asn1_indefinite_tag_end(asn1_data);
    asn1_indefinite_tag_end(asn1_data);
    d->date_type = 2;
  }

  return 1;
}

int get_db_tag(asn1_t *asn1_data, db_tag *tag) {
  bzero(tag, sizeof(db_tag));

  (void)asn1_sequence_start(asn1_data);
  asn1_get_mandatory_visible_string_field(asn1_data, &tag->db);
  get_object_id(asn1_data, &tag->tag);
  asn1_indefinite_tag_end(asn1_data);

  return 1;
}

int get_pdb_seq_id(asn1_t *asn1_data, pdb_seqid *seqid) {
  bzero(seqid, sizeof(pdb_seqid));

  (void)asn1_sequence_start(asn1_data);
  asn1_get_mandatory_visible_string_field(asn1_data, &seqid->mol);
  asn1_get_optional_integer_field(asn1_data, 0xA1, &seqid->chain);
  if (asn1_peek_byte(asn1_data) == 0xA2) {
    get_explicit_tag(asn1_data);
    get_date(asn1_data, &seqid->rel);
    asn1_indefinite_tag_end(asn1_data);
  }
  asn1_get_optional_visible_string_field(asn1_data, 0xA3, &seqid->chain_id);
  asn1_indefinite_tag_end(asn1_data);

  return 1;
}

int get_gi_import_id(asn1_t *asn1_data, gi_import_id *seqid) {
  bzero(seqid, sizeof(gi_import_id));

  (void)asn1_sequence_start(asn1_data);
  asn1_get_mandatory_integer_field(asn1_data, &seqid->id);
  asn1_get_optional_visible_string_field(asn1_data, 0xA1, &seqid->db);
  asn1_get_optional_visible_string_field(asn1_data, 0xA2, &seqid->release);
  asn1_indefinite_tag_end(asn1_data);

  return 1;
}

int get_text_seq_id(asn1_t *asn1_data, text_seqid *seqid) {
  bzero(seqid, sizeof(text_seqid));

  (void)asn1_sequence_start(asn1_data);
  asn1_get_optional_visible_string_field(asn1_data, 0xA0, &seqid->name);
  asn1_get_optional_visible_string_field(asn1_data, 0xA1, &seqid->acc);
  asn1_get_optional_visible_string_field(asn1_data, 0xA2, &seqid->rel);
  asn1_get_optional_integer_field(asn1_data, 0xA3, &seqid->ver);
  asn1_indefinite_tag_end(asn1_data);

  return 1;
}

int get_patent_seq_id(asn1_t *asn1_data, patent_seqid *seqid) {
  bzero(seqid, sizeof(patent_seqid));

  (void)asn1_sequence_start(asn1_data);
  asn1_get_mandatory_integer_field(asn1_data, &seqid->seqid);
  get_id_pat(asn1_data, seqid->cit);
  asn1_indefinite_tag_end(asn1_data);

  return 1;
}

int get_object_id(asn1_t *asn1_data, object_id *id) {
  bzero(id, sizeof(object_id));

  get_explicit_tag(asn1_data);
  uint8_t field_no = asn1_peek_byte(asn1_data);
  if (field_no == 0xA0) {
    asn1_get_mandatory_integer_field(asn1_data, &id->id.id);
    id->id_type = 1;
  } else {
    asn1_get_mandatory_visible_string_field(asn1_data, &id->id.str);
    id->id_type = 2;
  }
  asn1_indefinite_tag_end(asn1_data);

  return 1;
}

int get_id_pat(asn1_t *asn1_data, id_pat *id) {
  bzero(id, sizeof(id_pat));

  (void)asn1_sequence_start(asn1_data);
  asn1_get_optional_visible_string_field(asn1_data, 0xA0, &id->country);
  uint8_t field_no = asn1_peek_byte(asn1_data);
  if (field_no == 0xA1) {
    asn1_get_mandatory_visible_string_field(asn1_data, &id->doc_number.number);
    id->doc_number_type = 1;
  } else {
    asn1_get_mandatory_visible_string_field(asn1_data,
                                            &id->doc_number.app_number);
    id->doc_number_type = 2;
  }
  asn1_get_optional_visible_string_field(asn1_data, 0xA3, &id->doc_type);
  asn1_indefinite_tag_end(asn1_data);

  return 1;
}

int get_seq_id(asn1_t *asn1_data, seq_id *seqid) {
  uint8_t field_no = asn1_peek_byte(asn1_data);
  switch (field_no) {
  case 0xA0: {
    get_object_id(asn1_data, &seqid->id.obj_id);
    seqid->seq_id_type = seq_local;
    break;
  }
  case 0xA1: {
    get_explicit_tag(asn1_data);
    asn1_get_integer(asn1_data, &seqid->id.int_id);
    asn1_indefinite_tag_end(asn1_data);
    seqid->seq_id_type = seq_gibbsq;
    break;
  }
  case 0xA2: {
    get_explicit_tag(asn1_data);
    asn1_get_integer(asn1_data, &seqid->id.int_id);
    asn1_indefinite_tag_end(asn1_data);
    seqid->seq_id_type = seq_gibbmt;
    break;
  }
  case 0xA3: {
    get_explicit_tag(asn1_data);
    get_gi_import_id(asn1_data, &seqid->id.giim_id);
    seqid->seq_id_type = seq_giim;
    asn1_indefinite_tag_end(asn1_data);
    break;
  }
  case 0xA4: {
    get_explicit_tag(asn1_data);
    get_text_seq_id(asn1_data, &seqid->id.text_id);
    seqid->seq_id_type = seq_genbank;
    asn1_indefinite_tag_end(asn1_data);
    break;
  }
  case 0xA5: {
    get_explicit_tag(asn1_data);
    get_text_seq_id(asn1_data, &seqid->id.text_id);
    seqid->seq_id_type = seq_embl;
    asn1_indefinite_tag_end(asn1_data);
    break;
  }
  case 0xA6: {
    get_explicit_tag(asn1_data);
    get_text_seq_id(asn1_data, &seqid->id.text_id);
    seqid->seq_id_type = seq_pir;
    asn1_indefinite_tag_end(asn1_data);
    break;
  }
  case 0xA7: {
    get_explicit_tag(asn1_data);
    get_text_seq_id(asn1_data, &seqid->id.text_id);
    seqid->seq_id_type = seq_swissprot;
    asn1_indefinite_tag_end(asn1_data);
    break;
  }
  case 0xA8: {
    get_explicit_tag(asn1_data);
    get_patent_seq_id(asn1_data, &seqid->id.pat_id);
    seqid->seq_id_type = seq_patent;
    asn1_indefinite_tag_end(asn1_data);
    break;
  }
  case 0xA9: {
    get_explicit_tag(asn1_data);
    get_text_seq_id(asn1_data, &seqid->id.text_id);
    seqid->seq_id_type = seq_other;
    asn1_indefinite_tag_end(asn1_data);
    break;
  }
  case 0xAA: {
    get_explicit_tag(asn1_data);
    get_db_tag(asn1_data, &seqid->id.db_tag_id);
    seqid->seq_id_type = seq_general;
    asn1_indefinite_tag_end(asn1_data);
    break;
  }
  case 0xAB: {
    get_explicit_tag(asn1_data);
    asn1_get_integer(asn1_data, &seqid->id.int_id);
    asn1_indefinite_tag_end(asn1_data);
    seqid->seq_id_type = seq_gi;
    break;
  }
  case 0xAC: {
    get_explicit_tag(asn1_data);
    get_text_seq_id(asn1_data, &seqid->id.text_id);
    asn1_indefinite_tag_end(asn1_data);
    seqid->seq_id_type = seq_ddbj;
    break;
  }
  case 0xAD: {
    get_explicit_tag(asn1_data);
    get_text_seq_id(asn1_data, &seqid->id.text_id);
    asn1_indefinite_tag_end(asn1_data);
    seqid->seq_id_type = seq_prf;
    break;
  }
  case 0xAE: {
    get_explicit_tag(asn1_data);
    get_pdb_seq_id(asn1_data, &seqid->id.pdb_id);
    asn1_indefinite_tag_end(asn1_data);
    seqid->seq_id_type = seq_pdb;
    break;
  }
  case 0xAF: {
    get_explicit_tag(asn1_data);
    get_text_seq_id(asn1_data, &seqid->id.text_id);
    asn1_indefinite_tag_end(asn1_data);
    seqid->seq_id_type = seq_tpg;
    break;
  }
  case 0xB0: {
    get_explicit_tag(asn1_data);
    get_text_seq_id(asn1_data, &seqid->id.text_id);
    asn1_indefinite_tag_end(asn1_data);
    seqid->seq_id_type = seq_tpe;
    break;
  }
  case 0xB1: {
    get_explicit_tag(asn1_data);
    get_text_seq_id(asn1_data, &seqid->id.text_id);
    asn1_indefinite_tag_end(asn1_data);
    seqid->seq_id_type = seq_tpd;
    break;
  }
  case 0xB2: {
    get_explicit_tag(asn1_data);
    get_text_seq_id(asn1_data, &seqid->id.text_id);
    asn1_indefinite_tag_end(asn1_data);
    seqid->seq_id_type = seq_gpipe;
    break;
  }
  case 0xB3: {
    get_explicit_tag(asn1_data);
    get_text_seq_id(asn1_data, &seqid->id.text_id);
    asn1_indefinite_tag_end(asn1_data);
    seqid->seq_id_type = seq_named_annot_track;
    break;
  }
  default:
    seqid->seq_id_type = seq_none;
    break;
  }

  return seqid->seq_id_type;
}

int get_ints(asn1_t *asn1_data, int_vec *vec) {
  uint32_t value;
  asn1_sequence_start(asn1_data);
  for (;;) {
    if (asn1_peek_byte(asn1_data) == 0x02) {
      asn1_get_integer(asn1_data, &value);
      kv_push(integer, *vec, value);
    } else {
      break;
    }
  }
  asn1_indefinite_tag_end(asn1_data);

  return 1;
}

int get_memberships(asn1_t *asn1_data, int_vec *memberships) {
  if (asn1_peek_byte(asn1_data) != 0xA3) {
    return 0;
  }

  get_explicit_tag(asn1_data);
  get_ints(asn1_data, memberships);
  asn1_indefinite_tag_end(asn1_data);

  return 1;
}

int get_links(asn1_t *asn1_data, int_vec *links) {
  if (asn1_peek_byte(asn1_data) != 0xA4) {
    return 0;
  }

  get_explicit_tag(asn1_data);
  get_ints(asn1_data, links);
  asn1_indefinite_tag_end(asn1_data);

  return 1;
}


int get_other_info(asn1_t *asn1_data, int_vec *other_info) {
  if (asn1_peek_byte(asn1_data) != 0xA5) {
    return 0;
  }

  get_explicit_tag(asn1_data);
  get_ints(asn1_data, other_info);
  asn1_indefinite_tag_end(asn1_data);

  return 1;
}

int get_blast_defline(asn1_t *asn1_data, blast_defline *defline) {
  if (!asn1_sequence_start(asn1_data)) {
    return 0;
  }
  asn1_get_optional_visible_string_field(asn1_data, 0xA0, &defline->title);
  get_explicit_tag(asn1_data);
  asn1_sequence_start(asn1_data);
  int seq_id_type;
  for (;;) {
    seq_id id;
    seq_id_type = get_seq_id(asn1_data, &id);
    if (seq_id_type == seq_none) {
      break;
    }
    kv_push(seq_id, defline->seq_ids, id);
  }
  asn1_indefinite_tag_end(asn1_data);
  asn1_indefinite_tag_end(asn1_data);
  asn1_get_optional_integer_field(asn1_data, 0xA2, &defline->taxid);
  get_memberships(asn1_data, &defline->memberships);
  get_links(asn1_data, &defline->links);
  get_other_info(asn1_data, &defline->other_info);
  asn1_indefinite_tag_end(asn1_data);

  return 1;
}

blast_defline *init_blast_defline() {
  blast_defline *defline = alloc_memory(NULL, sizeof(blast_defline), 0, 1);
  kv_init(defline->seq_ids);
  kv_init(defline->memberships);
  kv_init(defline->links);
  kv_init(defline->other_info);

  return defline;
}

void reset_blast_defline(blast_defline *defline) {
  defline->title = NULL;
  kv_clear(defline->seq_ids);
  kv_clear(defline->memberships);
  kv_clear(defline->links);
  kv_clear(defline->other_info);
  defline->taxid = 0;

}

int get_blast_deflines(asn1_t *asn1_data, blast_deflines *deflines) {
  int retval;
  int i = 0;

  asn1_data->buf_len = 0;
  asn1_data->curr_pos = 0;

  if (!load_more_data(asn1_data)) {
    return 0;
  }
  asn1_sequence_start(asn1_data);
  do {
    blast_defline *defline = NULL;
    if (i == kv_size(*deflines)) {
      defline = init_blast_defline();
      kv_push(blast_defline *, *deflines, defline);
    } else {
      defline = kv_a(blast_defline *, *deflines, i);
    }
    reset_blast_defline(defline);
    retval = get_blast_defline(asn1_data, defline);
    if (retval == 0) {
      errx(1, "error parsing defline");
    }
    i++;
  } while (!asn1_is_end_of_sequence(asn1_data)) ;
  asn1_indefinite_tag_end(asn1_data);

  return i;
}
