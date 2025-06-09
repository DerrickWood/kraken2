#include <assert.h>
#include <err.h>
#include <fcntl.h>
#include <getopt.h>
#include <inttypes.h>
#include <libgen.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>

#include <sys/mman.h>
#include <netinet/in.h>

#include "blast_to_fasta.h"
#include "blast_defline.h"
#include "blast_utils.h"
#include "kvec.h"

/**
 * Mapping from masks to ASCII characters for ambiguous nucleotides.
 */
char mask2dna[] = {
        '?', // 0
        'A', // 1
        'C', // 2
        'M', // 3
        'G', // 4
        'R', // 5
        'S', // 6
        'V', // 7
        'T', // 8
        'W', // 9
        'Y', // 10
        'H', // 11
        'K', // 12
        'D', // 13
        'B', // 14
        'N', // 15 (inclusive N)
        'N'  // 16 (exclusive N)
};

const char *const four_mers[] = {
    "AAAA", "AAAC", "AAAG", "AAAT", "AACA", "AACC", "AACG", "AACT", "AAGA",
    "AAGC", "AAGG", "AAGT", "AATA", "AATC", "AATG", "AATT", "ACAA", "ACAC",
    "ACAG", "ACAT", "ACCA", "ACCC", "ACCG", "ACCT", "ACGA", "ACGC", "ACGG",
    "ACGT", "ACTA", "ACTC", "ACTG", "ACTT", "AGAA", "AGAC", "AGAG", "AGAT",
    "AGCA", "AGCC", "AGCG", "AGCT", "AGGA", "AGGC", "AGGG", "AGGT", "AGTA",
    "AGTC", "AGTG", "AGTT", "ATAA", "ATAC", "ATAG", "ATAT", "ATCA", "ATCC",
    "ATCG", "ATCT", "ATGA", "ATGC", "ATGG", "ATGT", "ATTA", "ATTC", "ATTG",
    "ATTT", "CAAA", "CAAC", "CAAG", "CAAT", "CACA", "CACC", "CACG", "CACT",
    "CAGA", "CAGC", "CAGG", "CAGT", "CATA", "CATC", "CATG", "CATT", "CCAA",
    "CCAC", "CCAG", "CCAT", "CCCA", "CCCC", "CCCG", "CCCT", "CCGA", "CCGC",
    "CCGG", "CCGT", "CCTA", "CCTC", "CCTG", "CCTT", "CGAA", "CGAC", "CGAG",
    "CGAT", "CGCA", "CGCC", "CGCG", "CGCT", "CGGA", "CGGC", "CGGG", "CGGT",
    "CGTA", "CGTC", "CGTG", "CGTT", "CTAA", "CTAC", "CTAG", "CTAT", "CTCA",
    "CTCC", "CTCG", "CTCT", "CTGA", "CTGC", "CTGG", "CTGT", "CTTA", "CTTC",
    "CTTG", "CTTT", "GAAA", "GAAC", "GAAG", "GAAT", "GACA", "GACC", "GACG",
    "GACT", "GAGA", "GAGC", "GAGG", "GAGT", "GATA", "GATC", "GATG", "GATT",
    "GCAA", "GCAC", "GCAG", "GCAT", "GCCA", "GCCC", "GCCG", "GCCT", "GCGA",
    "GCGC", "GCGG", "GCGT", "GCTA", "GCTC", "GCTG", "GCTT", "GGAA", "GGAC",
    "GGAG", "GGAT", "GGCA", "GGCC", "GGCG", "GGCT", "GGGA", "GGGC", "GGGG",
    "GGGT", "GGTA", "GGTC", "GGTG", "GGTT", "GTAA", "GTAC", "GTAG", "GTAT",
    "GTCA", "GTCC", "GTCG", "GTCT", "GTGA", "GTGC", "GTGG", "GTGT", "GTTA",
    "GTTC", "GTTG", "GTTT", "TAAA", "TAAC", "TAAG", "TAAT", "TACA", "TACC",
    "TACG", "TACT", "TAGA", "TAGC", "TAGG", "TAGT", "TATA", "TATC", "TATG",
    "TATT", "TCAA", "TCAC", "TCAG", "TCAT", "TCCA", "TCCC", "TCCG", "TCCT",
    "TCGA", "TCGC", "TCGG", "TCGT", "TCTA", "TCTC", "TCTG", "TCTT", "TGAA",
    "TGAC", "TGAG", "TGAT", "TGCA", "TGCC", "TGCG", "TGCT", "TGGA", "TGGC",
    "TGGG", "TGGT", "TGTA", "TGTC", "TGTG", "TGTT", "TTAA", "TTAC", "TTAG",
    "TTAT", "TTCA", "TTCC", "TTCG", "TTCT", "TTGA", "TTGC", "TTGG", "TTGT",
    "TTTA", "TTTC", "TTTG", "TTTT"};

string_t init_string() {
  string_t s;
  bzero(&s, sizeof(string_t));

  return s;
}

uint32_t string_length(string_t *s) { return s->len; }
uint32_t string_capacity(string_t *s) { return s->cap; }
void string_clear(string_t *s) { s->len = 0; }

uint8_t *string_data(string_t *s) { return s->string; }

uint8_t *to_c_string(string_t *s) {
  uint8_t *data = string_data(s);
  data[s->len] = '\0';

  return data;
}

void string_append_char(string_t *s, uint8_t c) {
  if (s->len == s->cap) {
    uint32_t new_cap = next_power_of_2(s->cap + 1);
    s->string = alloc_memory(s->string, 1, s->cap, new_cap);
    s->cap = new_cap;
  }
  s->string[s->len++] = c;
}

void string_append_int(string_t *s, uint32_t num) {
  int num_digits = (uint32_t)log10(num) + 1;
  int divisor = pow(10, num_digits - 1);

  while (num_digits-- > 0) {
    char digit = (num / divisor) + '0';
    num %= divisor;
    divisor /= 10;
    string_append_char(s, digit);
  }
}

void string_append_str(string_t *s, const char *str) {
  while (*str) {
    string_append_char(s, *str++);
  }
}

uint32_t string_copy(string_t *s1, string_t *s2, uint32_t s1_start,
                  uint32_t s2_start, uint32_t len) {
  if (s1_start > s1->len || s2_start > s2->len) {
    return 0;
  }

  len = len > s2->len ? s2->len : len;

  if (s1->cap < len) {
    uint32_t new_len = next_power_of_2(len + 1);
    s1->string = alloc_memory(s1->string, 1, s1->cap, new_len);
    s1->cap = new_len;
  }

  if (s1 == s2) {
    memmove(s1->string + s1_start, s2->string + s2_start, len);
  } else {
    memcpy(s1->string + s1_start, s2->string + s2_start, len);
  }

  s1->len = len + s1_start;

  return len;
}

uint32_t string_append(string_t *s1, string_t *s2, uint32_t s2_start,
                       uint32_t len) {
  return string_copy(s1, s2, s1->len, s2_start, len);
}

void string_reserve(string_t *s, uint32_t size) {
  s->string = alloc_memory(s->string, 1, s->cap, size + 1);

  if (size > s->cap) {
    s->cap = size + 1;
  }
}

void string_read_from_file(string_t *string, FILE *f, uint32_t offset,
                           uint32_t amt) {
  // Add error checking for offset being larger than string cap
  string_reserve(string, amt + offset);
  string->len = read_into_buffer(f, string->string + offset, 1, amt) + offset;
}

void free_string(string_t *s) {
  if (s->string != NULL)
    free(s->string);
  s->string = NULL;
  s->cap = 0;
  s->len = 0;
}

uint64_t ntoh_64(uint64_t val) {
  uint32_t a, b;

  a = (val >> 32) & 0xffffffff;
  b = val & 0xffffffff;

  val = ntohl(b);
  val = val << 32;
  val = val | ntohl(a);

  return val;
}

amb32_t u32_to_amb32(uint32_t n) { return *((amb32_t *)(&n)); }
amb64_t u64_to_amb64(uint64_t n) { return *((amb64_t *)(&n)); }

uint32_t read_int(FILE *f) {
  uint32_t value = 0;

  read_into_buffer(f, &value, sizeof(uint32_t), 1);
  return ntohl(value);
}

uint32_t read_long(FILE *f) {
  uint64_t value = 0;

  read_into_buffer(f, &value, sizeof(uint64_t), 1);
  return ntoh_64(value);
}

string_t read_string(FILE *f) {
  uint32_t len = 0;
  string_t string = init_string();

  len = read_int(f);
  string_reserve(&string, len + 1);

  read_into_buffer(f, string.string, 1, len);
  string.len = len;

  return string;
}

void *read_array(FILE *f, uint32_t element_size, uint32_t length) {
  void *array = alloc_memory(NULL, element_size, 0, length + 1);
  read_into_buffer(f, array, element_size, length);

  return array;
}

idx_t *init_idx_data(char *filename) {
  idx_t *idx = malloc(sizeof(idx_t));
  bzero(idx, sizeof(idx_t));

  idx->idx_file = open_file(filename, "r");

  return idx;
}

void read_idx_data(idx_t *idx_data) {
  idx_data->fmt_version = read_int(idx_data->idx_file);
  idx_data->db_seq_type = read_int(idx_data->idx_file);
  idx_data->volume = read_int(idx_data->idx_file);
  idx_data->title = read_string(idx_data->idx_file);
  idx_data->lmdb_file = read_string(idx_data->idx_file);
  idx_data->date = read_string(idx_data->idx_file);
  idx_data->num_oids = read_int(idx_data->idx_file);
  idx_data->vol_len = read_long(idx_data->idx_file);
  idx_data->max_seq_len = read_int(idx_data->idx_file);

  uint32_t array_bytes = 4 * (idx_data->num_oids + 1);

  idx_data->hdr_arr = read_array(idx_data->idx_file, 1, array_bytes);
  idx_data->seq_arr = read_array(idx_data->idx_file, 1, array_bytes);
  idx_data->amb_arr = read_array(idx_data->idx_file, 1, array_bytes);
}

void free_idx_data(idx_t *idx_data) {
  if (idx_data == NULL) {
    return;
  }

  if (idx_data->idx_file != NULL) {
    fclose(idx_data->idx_file);
  }

  free_string(&idx_data->title);
  free_string(&idx_data->lmdb_file);
  free_string(&idx_data->date);

  if (idx_data->hdr_arr != NULL) {
    free(idx_data->hdr_arr);
  }
  if (idx_data->seq_arr != NULL) {
    free(idx_data->seq_arr);
  }
  if (idx_data->amb_arr != NULL) {
    free(idx_data->amb_arr);
  }

  free(idx_data);
}

hdr_t *init_hdr_data(char *filename, uint32_t *hdr_offsets, uint32_t num_oids) {
  hdr_t *hdr_data = malloc(sizeof(hdr_t));
  bzero(hdr_data, sizeof(hdr_t));

  /* hdr_data->hdr_file = open_file(filename, "r"); */
  hdr_data->fasta_hdr = init_string();
  kv_init(hdr_data->deflines);
  for (uint32_t i = 1; i < num_oids + 1; i++) {
    hdr_offsets[i - 1] = ntohl(hdr_offsets[i]) - ntohl(hdr_offsets[i - 1]);
  }
  hdr_data->asn1 = init_asn1(filename, hdr_offsets, num_oids);
  hdr_offsets[num_oids] = 0;

  return hdr_data;
}


void free_hdr_data(hdr_t *hdr) {
  if (hdr == NULL) {
    return;
  }

  if (hdr->hdr_file != NULL) {
    fclose(hdr->hdr_file);
  }
  kv_destroy(hdr->deflines);
  free_string(&hdr->fasta_hdr);
  free_asn1(hdr->asn1);

  free(hdr);
  hdr = NULL;
}

void db_tag_to_string(string_t *s, db_tag tag) {
  string_append_str(s, (const char *)tag.db);
  string_append_char(s, '_');
  if (tag.tag.id_type == 1) {
    string_append_int(s, tag.tag.id.id);
  } else {
    visible_string str = tag.tag.id.str;
    string_append_str(s, (const char *)str);
  }
}

void text_seq_id_to_string(string_t *s, text_seqid seqid) {
  visible_string acc = seqid.acc;
  string_append_str(s, (const char *)acc);

  if (seqid.ver > 0) {
    string_append_char(s, '.');
    string_append_int(s, seqid.ver);
  }
}

void pdb_seq_id_to_string(string_t *s, pdb_seqid seqid) {
  visible_string str = seqid.mol;
  string_append_str(s, (const char *)str);

  str = seqid.chain_id;
  if (str) {
    string_append_char(s, '_');
    string_append_str(s, (const char *)str);
  } else if (seqid.chain != 0) {
    string_append_char(s, '_');
    string_append_char(s, (char)seqid.chain);
  }
}


void defline_to_header(string_t *s, blast_deflines deflines, int include_taxid, int curr_defline) {
  string_append_char(s, '>');
  blast_defline *defline = kv_A(deflines, curr_defline);

  if (include_taxid && defline->taxid) {
    string_append_str(s, "kraken:taxid|");
    string_append_int(s, defline->taxid);
    string_append_char(s, '|');
  }

  for (int i = 0; i < kv_size(defline->seq_ids); i++) {
    switch (kv_A(defline->seq_ids, i).seq_id_type) {
    case seq_gi:
      break;
    case seq_pdb:
      pdb_seq_id_to_string(s, kv_A(defline->seq_ids, i).id.pdb_id);
      break;
    case seq_general:
      db_tag_to_string(s, kv_A(defline->seq_ids, i).id.db_tag_id);
      break;
    default:
      text_seq_id_to_string(s, kv_A(defline->seq_ids, i).id.text_id);
    }
  }

  string_append_char(s, ' ');
  visible_string title = defline->title;
  string_append_str(s, (const char *)title);
}

void deflines_to_header(string_t *s, blast_deflines deflines, int include_taxid, int num_deflines) {
  for (int i = 0; i < num_deflines; i++) {
    if (string_length(s) > 0) {
      string_append_char(s, ' ');
    }
    defline_to_header(s, deflines, include_taxid, i);
  }
}

int get_deflines(hdr_t *hdr_data) {
  int num_deflines = get_blast_deflines(hdr_data->asn1, &hdr_data->deflines);
  if (num_deflines == 0) {
    return 0;
  }

  return num_deflines;
}

seq_t *init_seq_data(char *filename, uint32_t max_seq_len, uint32_t max_buf_size) {
  char first_byte;

  seq_t *seq_data = malloc(sizeof(seq_t));
  bzero(seq_data, sizeof(seq_t));

  seq_data->seq_file = open_file(filename, "r");
  // discard the first byte of this file
  fread(&first_byte, 1, 1, seq_data->seq_file);
  seq_data->seq = init_string();
  string_reserve(&seq_data->seq, max_seq_len + 1);

  seq_data->buffer = init_string();
  string_reserve(&seq_data->buffer, max_buf_size + 1);

  return seq_data;
}


uint32_t max_block_size(uint32_t *offsets, uint32_t length) {
  uint32_t max_length = 0;

  for (uint32_t i = 1; i < length; i++) {
    uint32_t delta = (ntohl(offsets[i]) - ntohl(offsets[i - 1]));
    if (delta > max_length)
      max_length = delta;
  }

  return max_length;
}

int has_ambiguous_data(uint32_t block_end, uint32_t amb_start) {
  return block_end != amb_start;
}

uint32_t get_nucleotide_length(uint8_t *data, size_t data_length) {
  // if the number of nucleotides takes up less than a byte
  // the last byte encodes both the nucleotides and the count
  uint8_t remainder = ntohl(data[data_length]);

  return (data_length - 1) * 4 + remainder;
}

void read_seq_block(seq_t *seq, uint32_t block_length) {
  string_clear(&seq->buffer);
  read_into_buffer(seq->seq_file, string_data(&seq->buffer), 1, block_length);
  seq->buffer.len = block_length;
}

uint32_t read_ambiguous_data(uint8_t *input, amb_t **output, int output_length, uint32_t *format) {
  uint32_t amb_data_len = 0;
  memcpy(&amb_data_len, input, sizeof(uint32_t));
  amb_data_len = ntohl(amb_data_len);
  *format = (amb_data_len & 0x80000000) == 0x80000000;
  amb_data_len = amb_data_len & 0x7fffffff;
  *output = alloc_memory(*output, sizeof(amb_t), output_length, amb_data_len);

  input += sizeof(uint32_t);
  if (*format == 1) {
    amb_data_len /= 2;
    uint64_t segment;
    amb64_t *out = *((amb64_t **)output);
    for (size_t i = 0; i < amb_data_len; i++) {
      memcpy(&segment, input, sizeof(uint64_t));
      out[i] = u64_to_amb64(ntoh_64(segment));
      input += sizeof(uint64_t);
    }
  } else {
    uint32_t segment;
    amb32_t *out = *((amb32_t **)output);
    for (size_t i = 0; i < amb_data_len; i++) {
      memcpy(&segment, input, sizeof(uint32_t));
      out[i] = u32_to_amb32(ntohl(segment));
      input += sizeof(uint32_t);
    }
  }

  return amb_data_len;
}

void reconstruct_sequence(seq_t *seq, uint32_t unamb_data_len, uint32_t amb_data_len, int format) {
  uint8_t *buffer = string_data(&seq->buffer);
  uint8_t final_byte = buffer[--unamb_data_len];
  uint32_t nucs_in_final_byte = final_byte & 0x3;
  uint32_t nuc_len = unamb_data_len * 4 + nucs_in_final_byte;

  int k = 3;
  uint32_t j = 0;
  const uint8_t bases[] = {'A', 'C', 'G', 'T'};
  uint8_t *seq_buf = string_data(&seq->seq);
  // decode every nucs in every byte except the last
  for (uint32_t i = 0; i < unamb_data_len; i++) {
    const char *kmer = four_mers[buffer[i]];
    seq_buf[j++] = *kmer++;
    seq_buf[j++] = *kmer++;
    seq_buf[j++] = *kmer++;
    seq_buf[j++] = *kmer++;
  }

  if (final_byte != 0) {
    while (nucs_in_final_byte-- > 0) {
      seq_buf[j] = bases[(final_byte >> (k-- << 1)) & 0x3];
      j++;
    }
  }

  if (format == 1) {
    amb64_t *amb_data = (amb64_t *)seq->amb_data;
    for (uint32_t i = 0; i < amb_data_len; i++) {
      amb64_t a = amb_data[i];
      // if (a.value != 0 && a.offset != 0 && a.value != 0) {
      memset(seq->seq.string + a.offset, mask2dna[a.value], a.length + 1);
      // }
    }
  } else {
    amb32_t *amb_data = (amb32_t *)seq->amb_data;
    for (uint32_t i = 0; i < amb_data_len; i++) {
      amb32_t a = amb_data[i];
      // if (a.value != 0 && a.offset != 0 && a.value != 0) {
      memset(seq->seq.string + a.offset, mask2dna[a.value], a.length + 1);
      // }
    }
  }

  assert(j == nuc_len);
  seq_buf[nuc_len] = '\0';
  seq->seq.len = nuc_len;
}

int next_sequence(seq_t *seq, idx_t *idx_data) {
  uint32_t i = seq->curr_pos;
  uint32_t format = 0;
  uint32_t amb_data_len = 0;
  uint32_t seq_start = ntohl(idx_data->seq_arr[i]);
  uint32_t seq_end = ntohl(idx_data->seq_arr[i + 1]);
  uint32_t amb_start = ntohl(idx_data->amb_arr[i]);

  read_seq_block(seq, seq_end - seq_start);
  if (has_ambiguous_data(seq_end, amb_start)) {
    uint8_t *buffer = string_data(&seq->buffer);
    amb_data_len =
      read_ambiguous_data(buffer + (amb_start - seq_start),
                             &seq->amb_data, seq->amb_data_cap, &format);
    if (amb_data_len > seq->amb_data_cap) {
      seq->amb_data_cap = amb_data_len;
    }
  }

  uint32_t unamb_data_len = amb_start - seq_start;
  reconstruct_sequence(seq, unamb_data_len, amb_data_len, format);
  seq->curr_pos = seq->curr_pos + 1;

  return 0;
}

void write_sequence(string_t *seq, int seq_width, FILE *f) {
  uint8_t *raw_seq = string_data(seq);
  int full_width_lines = seq->len / seq_width;
  int remainder = seq->len % seq_width;

  for (int i = 0; i < full_width_lines; i++) {
    fwrite(raw_seq, 1, seq_width, f);
    fputc('\n', f);
    raw_seq += seq_width;
  }

  if (remainder) {
    fwrite(raw_seq, 1, remainder, f);
    fputc('\n', f);
  }
}

void free_seq_data(seq_t *seq_data) {
  if (seq_data == NULL) {
    return;
  }

  if (seq_data->seq_file != NULL) {
    fclose(seq_data->seq_file);
  }

  free_string(&seq_data->buffer);
  free_string(&seq_data->seq);

  if (seq_data->amb_data != NULL) {
    free(seq_data->amb_data);
  }

  free(seq_data);
}

void usage(const char *prog) {
  fprintf(stderr, "%s [-hst] [-o out_file] [-w width] blast_volume\n", prog);
  exit(EXIT_SUCCESS);
}

void help(const char *prog) {
  fprintf(stderr, "%s [-hst] [-o out_file] [-w width] blast_volume\n", prog);
  fprintf(stderr, "  blast_volume: the base name of the blast volume e.g. core_nt.00\n");
  fprintf(stderr, "  -h: print this help message and exit\n");
  fprintf(stderr, "  -o: The filename that the output gets saved to (default: "
                  "<blast_volume>.fna)\n");
  fprintf(stderr, "  -s: BLAST merges the headers of FASTA entries with "
                  "identical sequences, this option outputs a complete FASTA record"
                  "for every such header\n");
  fprintf(stderr, "  -t: Prepend the tax ID to the FASTA header with format "
                  "kraken:taxid|12345|\n");
  fprintf(stderr, "  -w: The width of the FASTA sequences (default: 80)\n");

  exit(EXIT_SUCCESS);
}

int main(int argc, char **argv) {
  char *volume;
  int option;
  int seq_width = 80;
  int split_hdr = 0;
  int include_taxid = 0;
  int out_filename_seen = 0;
  char *out_filename = NULL;
  const char *prog = basename(argv[0]);

  while ((option = getopt_long(argc, argv, "ho:stw:", NULL, NULL)) != -1) {
    switch (option) {
    case 'h':
      help(argv[0]);
    case 'o': {
      out_filename = optarg;
      out_filename_seen = 1;
      break;
    }
    case 's':
      split_hdr = 1;
      break;
    case 't':
      include_taxid = 1;
      break;
    case 'w':
      seq_width = strtol(optarg, NULL, 10);
      if (seq_width < 1) {
        fprintf(stderr, "-w has to be at least 1\n");
        exit(EXIT_FAILURE);
      }
      break;
    case '?':
    default:
      usage(argv[0]);
    }
  }
  argc -= optind;
  argv += optind;

  if (argc != 1) {
    fprintf(stderr, "Missing basename of BLAST volume\n");
    usage(prog);
  }

  volume = argv[0];
  int ext_len = 4;
  int volume_len = strlen(volume);
  char *idx_filename = alloc_memory(NULL, 1, 0, volume_len + ext_len + 1);
  char *hdr_filename = alloc_memory(NULL, 1, 0, volume_len + ext_len + 1);
  char *seq_filename = alloc_memory(NULL, 1, 0, volume_len + ext_len + 1);

  sprintf(idx_filename, "%s.nin", volume);
  sprintf(hdr_filename, "%s.nhr", volume);
  sprintf(seq_filename, "%s.nsq", volume);

  if (!out_filename_seen) {
    out_filename = alloc_memory(NULL, 1, 0, volume_len + ext_len + 1);
    sprintf(out_filename, "%s.fna", volume);
  }

  idx_t *idx_data = init_idx_data(idx_filename);
  read_idx_data(idx_data);
  hdr_t *hdr_data = init_hdr_data(hdr_filename, idx_data->hdr_arr, idx_data->num_oids);
  uint32_t max_seq_block_size =
      max_block_size(idx_data->seq_arr, idx_data->num_oids);
  seq_t *seq_data = init_seq_data(seq_filename, idx_data->max_seq_len, max_seq_block_size);

  fprintf(stderr, "The version of this database is %X\n", idx_data->fmt_version);
  fprintf(stderr, "Sequence type is: %s\n", idx_data->db_seq_type == 0 ? "nucleotide" : "protein");
  fprintf(stderr, "Volume number: %d\n", idx_data->volume);
  fprintf(stderr, "Title of volume: %s\n", to_c_string(&idx_data->title));
  fprintf(stderr, "Date created: %s\n", to_c_string(&idx_data->date));
  fprintf(stderr, "Number of OIDs: %d\n", idx_data->num_oids);
  fprintf(stderr, "Maximum sequence length: %d\n", idx_data->max_seq_len);

  FILE *out_file = open_file(out_filename, "w");

  for (uint32_t i = 0; i < idx_data->num_oids + 1; i++) {
    int num_headers = get_deflines(hdr_data);
    if (num_headers == 0) {
      break;
    }
    next_sequence(seq_data, idx_data);
    if (split_hdr) {
      for (int i = 0; i < num_headers; i++) {
        string_clear(&hdr_data->fasta_hdr);
        defline_to_header(&hdr_data->fasta_hdr, hdr_data->deflines, include_taxid, i);
        fprintf(out_file, "%s\n", to_c_string(&hdr_data->fasta_hdr));
        write_sequence(&seq_data->seq, seq_width, out_file);
      }
    } else {
      string_clear(&hdr_data->fasta_hdr);
      deflines_to_header(&hdr_data->fasta_hdr, hdr_data->deflines,
                         include_taxid, num_headers);
      fprintf(out_file, "%s\n", to_c_string(&hdr_data->fasta_hdr));
      write_sequence(&seq_data->seq, seq_width, out_file);
    }
  }

  fclose(out_file);
  free_idx_data(idx_data);
  free_hdr_data(hdr_data);
  free_seq_data(seq_data);

  free(idx_filename);
  free(hdr_filename);
  free(seq_filename);
  if (!out_filename_seen) {
    free(out_filename);
  }

  return 0;
}
