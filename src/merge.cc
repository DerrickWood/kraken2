#include "kraken2_data.h"
#include "reports.h"
#include "taxonomy.h"
#include <cstdlib>
#include <cstring>
#include <err.h>
#include <string.h>
#include <unistd.h>
#include <vector>
#include <stdint.h>

enum merge_errors {
        no_error = 0,
        no_output,
        no_taxon,
        no_input,
        one_input,
        too_much_input,
        invalid_flag,
};

typedef int taxid_t;

struct taxid_and_count {
        int taxid;
        int count;
};

static const int AMBIGUOUS_TAXID = INT32_MAX;
static const int MATE_PAIR_BORDER_TAXON = INT32_MAX - 1;
static const int READING_FRAME_BORDER_TAXON = INT32_MAX - 2;
static const int VALID_TAXID_THRESHOLD = INT32_MAX - 3;

#define __FJ__FALLTHROUGH (void)0


// credit: Paul Khuong
uint64_t encode_ten_thousands(uint64_t hi, uint64_t lo) {
        uint64_t merged = hi | (lo << 32);
        uint64_t top = ((merged * 10486ULL) >> 20) & ((0x7FULL << 32) | 0x7FULL);
        uint64_t bot = merged - 100ULL * top;
        uint64_t hundreds;
        uint64_t tens;
        hundreds = (bot << 16) + top;
        tens = (hundreds * 103ULL) >> 10;
        tens &= (0xFULL << 48) | (0xFULL << 32) | (0xFULL << 16) | 0xFULL;
        tens += (hundreds - 10ULL * tens) << 8;

        return tens;
}

char *to_string_khuong(uint64_t x, char *out) {
        if (x < 10) {}
        uint64_t top = x / 100000000;
        uint64_t bottom = x % 100000000;
        uint64_t first =
                0x3030303030303030 + encode_ten_thousands(top / 10000, top % 10000);
        memcpy(out, &first, sizeof(first));
        uint64_t second =
                0x3030303030303030 + encode_ten_thousands(bottom / 10000, bottom % 10000);
        memcpy(out + 8, &second, sizeof(second));

        while (*out == '0')
          out += 1;

        return out;
}

// credit: https://gist.github.com/niXman/5c0e53ad0dc98e66399658915747828e#file-strtoint-cpp-L90
std::uint64_t nixmans_atou64_shift(const char *ptr, std::size_t len) {
#   define __FJ__PER_CHAR_EXPR(n) ((res << 1) + (res << 3) + (str[len - n] - '0'))

        const auto *str = reinterpret_cast<const std::uint8_t *>(ptr);
        std::uint64_t res = 0;
        switch ( len ) {
        case 20: res = __FJ__PER_CHAR_EXPR(20); __FJ__FALLTHROUGH;
        case 19: res = __FJ__PER_CHAR_EXPR(19); __FJ__FALLTHROUGH;
        case 18: res = __FJ__PER_CHAR_EXPR(18); __FJ__FALLTHROUGH;
        case 17: res = __FJ__PER_CHAR_EXPR(17); __FJ__FALLTHROUGH;
        case 16: res = __FJ__PER_CHAR_EXPR(16); __FJ__FALLTHROUGH;
        case 15: res = __FJ__PER_CHAR_EXPR(15); __FJ__FALLTHROUGH;
        case 14: res = __FJ__PER_CHAR_EXPR(14); __FJ__FALLTHROUGH;
        case 13: res = __FJ__PER_CHAR_EXPR(13); __FJ__FALLTHROUGH;
        case 12: res = __FJ__PER_CHAR_EXPR(12); __FJ__FALLTHROUGH;
        case 11: res = __FJ__PER_CHAR_EXPR(11); __FJ__FALLTHROUGH;
        case 10: res = __FJ__PER_CHAR_EXPR(10); __FJ__FALLTHROUGH;
        case 9 : res = __FJ__PER_CHAR_EXPR( 9); __FJ__FALLTHROUGH;
        case 8 : res = __FJ__PER_CHAR_EXPR( 8); __FJ__FALLTHROUGH;
        case 7 : res = __FJ__PER_CHAR_EXPR( 7); __FJ__FALLTHROUGH;
        case 6 : res = __FJ__PER_CHAR_EXPR( 6); __FJ__FALLTHROUGH;
        case 5 : res = __FJ__PER_CHAR_EXPR( 5); __FJ__FALLTHROUGH;
        case 4 : res = __FJ__PER_CHAR_EXPR( 4); __FJ__FALLTHROUGH;
        case 3 : res = __FJ__PER_CHAR_EXPR( 3); __FJ__FALLTHROUGH;
        case 2 : res = __FJ__PER_CHAR_EXPR( 2); __FJ__FALLTHROUGH;
        case 1 : res = __FJ__PER_CHAR_EXPR( 1);
        }

#   undef __FJ__PER_CHAR_EXPR

        return res;
}

// https://lemire.me/blog/2021/11/18/converting-integers-to-fix-digit-representations-quickly/
char *int_to_string(uint64_t x, char *out) {
        static const char table[200] = {
                0x30, 0x30, 0x30, 0x31, 0x30, 0x32, 0x30, 0x33, 0x30, 0x34, 0x30, 0x35,
                0x30, 0x36, 0x30, 0x37, 0x30, 0x38, 0x30, 0x39, 0x31, 0x30, 0x31, 0x31,
                0x31, 0x32, 0x31, 0x33, 0x31, 0x34, 0x31, 0x35, 0x31, 0x36, 0x31, 0x37,
                0x31, 0x38, 0x31, 0x39, 0x32, 0x30, 0x32, 0x31, 0x32, 0x32, 0x32, 0x33,
                0x32, 0x34, 0x32, 0x35, 0x32, 0x36, 0x32, 0x37, 0x32, 0x38, 0x32, 0x39,
                0x33, 0x30, 0x33, 0x31, 0x33, 0x32, 0x33, 0x33, 0x33, 0x34, 0x33, 0x35,
                0x33, 0x36, 0x33, 0x37, 0x33, 0x38, 0x33, 0x39, 0x34, 0x30, 0x34, 0x31,
                0x34, 0x32, 0x34, 0x33, 0x34, 0x34, 0x34, 0x35, 0x34, 0x36, 0x34, 0x37,
                0x34, 0x38, 0x34, 0x39, 0x35, 0x30, 0x35, 0x31, 0x35, 0x32, 0x35, 0x33,
                0x35, 0x34, 0x35, 0x35, 0x35, 0x36, 0x35, 0x37, 0x35, 0x38, 0x35, 0x39,
                0x36, 0x30, 0x36, 0x31, 0x36, 0x32, 0x36, 0x33, 0x36, 0x34, 0x36, 0x35,
                0x36, 0x36, 0x36, 0x37, 0x36, 0x38, 0x36, 0x39, 0x37, 0x30, 0x37, 0x31,
                0x37, 0x32, 0x37, 0x33, 0x37, 0x34, 0x37, 0x35, 0x37, 0x36, 0x37, 0x37,
                0x37, 0x38, 0x37, 0x39, 0x38, 0x30, 0x38, 0x31, 0x38, 0x32, 0x38, 0x33,
                0x38, 0x34, 0x38, 0x35, 0x38, 0x36, 0x38, 0x37, 0x38, 0x38, 0x38, 0x39,
                0x39, 0x30, 0x39, 0x31, 0x39, 0x32, 0x39, 0x33, 0x39, 0x34, 0x39, 0x35,
                0x39, 0x36, 0x39, 0x37, 0x39, 0x38, 0x39, 0x39,
        };

        if (x > VALID_TAXID_THRESHOLD) {
                switch (x) {
                case AMBIGUOUS_TAXID:
                        out[0] = 'A';
                        break;
                case MATE_PAIR_BORDER_TAXON:
                        out[0] = '|';
                        break;
                case READING_FRAME_BORDER_TAXON:
                        out[0] = '-';
                        break;
                default:
                        break;
                }
                out[1] = '\0';

                return out;
        }

        uint64_t top = x / 100000000;
        uint64_t bottom = x % 100000000;
        //
        uint64_t toptop = top / 10000;
        uint64_t topbottom = top % 10000;
        uint64_t bottomtop = bottom / 10000;
        uint64_t bottombottom = bottom % 10000;
        //
        uint64_t toptoptop = toptop / 100;
        uint64_t toptopbottom = toptop % 100;

        uint64_t topbottomtop = topbottom / 100;
        uint64_t topbottombottom = topbottom % 100;

        uint64_t bottomtoptop = bottomtop / 100;
        uint64_t bottomtopbottom = bottomtop % 100;

        uint64_t bottombottomtop = bottombottom / 100;
        uint64_t bottombottombottom = bottombottom % 100;
        //
        memcpy(out, &table[2 * toptoptop], 2);
        memcpy(out + 2, &table[2 * toptopbottom], 2);
        memcpy(out + 4, &table[2 * topbottomtop], 2);
        memcpy(out + 6, &table[2 * topbottombottom], 2);
        memcpy(out + 8, &table[2 * bottomtoptop], 2);
        memcpy(out + 10, &table[2 * bottomtopbottom], 2);
        memcpy(out + 12, &table[2 * bottombottomtop], 2);
        memcpy(out + 14, &table[2 * bottombottombottom], 2);

        size_t i = 0;
        while (out[i] == '0' && i < 15)
                i += 1;
        out[16] = '\0';

        return &out[i];
}

void write_hit_list(vector<taxid_and_count> &hit_list, const char *border, FILE *f) {
        if (hit_list.empty()) {
                return;
        }

        char itoa_buffer[17] = { 0 };
        int previous_taxid = hit_list[0].taxid;
        int previous_count = hit_list[0].count;

        for (size_t i = 1; i < hit_list.size(); i++) {
                int taxid = hit_list[i].taxid;
                int count = hit_list[i].count;

                if (taxid == 0 && count == 0) {
                        fputs(int_to_string(previous_taxid, itoa_buffer), f);
                        fputc(':', f);
                        fputs(int_to_string(previous_count, itoa_buffer), f);
                        previous_taxid = taxid;
                        previous_count = count;
                }

                if (previous_taxid == taxid) {
                        previous_count += count;
                } else if (previous_count == 0 && count > 0) {
                        previous_count = count;
                        previous_taxid = taxid;
                } else {
                        fputs(int_to_string(previous_taxid, itoa_buffer), f);
                        fputc(':', f);
                        fputs(int_to_string(previous_count, itoa_buffer), f);
                        fputc(' ', f);
                        previous_taxid = taxid;
                        previous_count = count;
                }
        }

        fputs(int_to_string(previous_taxid, itoa_buffer), f);
        fputc(':', f);
        fputs(int_to_string(previous_count, itoa_buffer), f);
}

int resolve_tree(kraken2::Taxonomy &taxonomy, kraken2::taxon_counts_t &hit_counts,
                 size_t total_minimizers, float confidence_threshold)
{
        int max_taxon = 0;
        uint32_t max_score = 0;
        uint32_t required_score = ceil(confidence_threshold * total_minimizers);

        // Sum each taxon's LTR path, find taxon with highest LTR score
        for (auto &kv_pair : hit_counts) {
                taxid_t taxon = taxonomy.GetInternalID(kv_pair.first);
                uint32_t score = 0;

                for (auto &kv_pair2 : hit_counts) {
                        taxid_t taxon2 = taxonomy.GetInternalID(kv_pair2.first);

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
        max_score = hit_counts[taxonomy.nodes()[max_taxon].external_id];
        // We probably have a call w/o required support (unless LCA resolved tie)
        while (max_taxon && max_score < required_score) {
                max_score = 0;
                for (auto &kv_pair : hit_counts) {
                        taxid_t taxon = taxonomy.GetInternalID(kv_pair.first);
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

        return taxonomy.nodes()[max_taxon].external_id;
}

void parse_hit_list(char *string, int len,
                    std::vector<taxid_and_count> &counts) {
        char *end_p = &string[len];
        size_t n_items = 0;
        size_t i = 0;

        for (size_t i = 0; i < len; i++) {
                n_items += string[i] == ':';
        }

        counts.clear();
        counts.resize(n_items);

        char *str_taxid;
        char *str_count;

        int count = 0;
        int taxid = 0;

        while (true) {
                size_t len1, len2;
                str_taxid = strsep(&string, ":");
                len1 = string - str_taxid - 1;
                str_count = strsep(&string, " ");
                len2 = (string == NULL ? end_p : string) - str_count - 1;

                if (!str_taxid || !str_count) {
                        break;
                }
                // count = strtol(str_count, NULL, 10);
                // taxid = strtol(str_taxid, NULL, 10);
                if (str_taxid[0] == 'A') {
                        taxid = AMBIGUOUS_TAXID;
                        count = nixmans_atou64_shift(str_count, len2);
                } else if (str_taxid[0] == '|') {
                        taxid = MATE_PAIR_BORDER_TAXON;
                        count = MATE_PAIR_BORDER_TAXON;
                } else if (str_taxid[0] == '-') {
                        taxid = READING_FRAME_BORDER_TAXON;
                        count = READING_FRAME_BORDER_TAXON;
                } else {
                        taxid = nixmans_atou64_shift(str_taxid, len1);
                        count = nixmans_atou64_shift(str_count, len2);

                }
                counts[i++] = {taxid, count};
        }
}

int get_lca(kraken2::Taxonomy &taxonomy, uint64_t taxid1, uint64_t taxid2) {
        taxid1 = taxonomy.GetInternalID(taxid1);
        taxid2 = taxonomy.GetInternalID(taxid2);

        uint64_t lca = taxonomy.LowestCommonAncestor(taxid1, taxid2);

        return (int)taxonomy.nodes()[lca].external_id;
}

size_t merge_hit_lists(kraken2::Taxonomy &taxonomy,
                       kraken2::taxon_counts_t &hit_counts,
                       std::vector<taxid_and_count> &hit_list1,
                       std::vector<taxid_and_count> &hit_list2,
                       std::vector<taxid_and_count> &merged_hit_list) {
        size_t i1 = 0;
        size_t i2 = 0;
        size_t total_minimizers = 0;
        // ambiguous taxons or separator identifiers do not count as
        // valid tax id. We should not add them to hit counts to
        // avoid any potential issues when resolving final taxid.
        bool add_to_hit_counts = true;

        while (true) {
                taxid_and_count &tc1 = hit_list1[i1];
                taxid_and_count &tc2 = hit_list2[i2];
                int final_taxid = 0;
                int final_count = 0;

                // We have either encountered a ambigous taxid
                // or a pair or translated search delimiter.
                // Do not bother running the LCA on those,
                // they will be used later when outputting the
                // hit list.
                if (tc1.taxid > VALID_TAXID_THRESHOLD) {
                        final_taxid = tc1.taxid;
                        add_to_hit_counts = false;
                } else {
                        final_taxid = get_lca(taxonomy, tc1.taxid, tc2.taxid);
                        add_to_hit_counts = true;
                }

                if (tc1.count < tc2.count) {
                        tc2.count -= tc1.count;
                        i1 += 1;
                        final_count = tc1.count;
                        i2 += tc2.count == 0;
                } else if (tc2.count < tc1.count) {
                        tc1.count -= tc2.count;
                        i2 += 1;
                        final_count = tc2.count;
                        i1 += tc1.count == 0;
                } else {
                        i1 += 1;
                        i2 += 1;
                        final_count = tc1.count;
                }

                total_minimizers += final_count;
                merged_hit_list.push_back({final_taxid, final_count});
                if (add_to_hit_counts) {
                        hit_counts[final_taxid] += final_count;
                }
                if (i1 >= hit_list1.size()) {
                        break;
                }
        }

        return total_minimizers;
}

void get_fields(char *line, const char *delim, char *fields[], int fields_len) {
        char **fp;

        for (fp = fields; (*fp = strsep(&line, delim)) != NULL;) {
                if (**fp != '\0' && ++fp >= &fields[fields_len]) {
                        break;
                }
        }
}

std::tuple<size_t, size_t>
merge_classification_output(kraken2::Taxonomy &taxonomy, FILE *in1, FILE *in2,
                            FILE *out, float confidence_threshold,
                            kraken2::taxon_counters_t *counters,
                            FILE *classified_headers) {
        char *line1 = NULL;
        char *line2 = NULL;

        size_t line1_cap = 0;
        size_t line2_cap = 0;

        ssize_t line1_len = 0;
        ssize_t line2_len = 0;

        char *fields1[5];
        char *fields2[5];

        const char *status;
        // const char *header;
        const char *taxid;
        const char *hit_list;
        // const char *seq_len;

        enum {
                status_field = 0,
                header_field,
                taxid_field,
                len_field,
                hit_list_field,
        };

        std::vector<taxid_and_count> hit_list1;
        std::vector<taxid_and_count> hit_list2;
        std::vector<taxid_and_count> merged_hit_list;
        kraken2::taxon_counts_t hit_counts;
        char itoa_buf[17];

        size_t total_sequences = 0;
        size_t total_unclassified = 0;

        while (true) {
                line1_len = getline(&line1, &line1_cap, in1);
                line2_len = getline(&line2, &line2_cap, in2);

                if (line1_len == -1) {
                        break;
                }

                get_fields(line1, "\t", fields1, 5);
                get_fields(line2, "\t", fields2, 5);

                if (fields1[status_field][0] == 'C' && fields2[status_field][0] == 'C') {
                        parse_hit_list(fields1[hit_list_field], strlen(fields1[hit_list_field]), hit_list1);
                        parse_hit_list(fields2[hit_list_field], strlen(fields2[hit_list_field]), hit_list2);

                        size_t total_minimizers =
                                merge_hit_lists(taxonomy, hit_counts, hit_list1, hit_list2, merged_hit_list);
                        int res = resolve_tree(taxonomy, hit_counts, total_minimizers, confidence_threshold);
                        taxid = int_to_string(res, itoa_buf);
                        status = "C";
                } else if (*fields1[status_field] == 'C' &&
                           *fields2[status_field] == 'U') {
                        status = "C";
                        taxid = fields1[taxid_field];
                        hit_list = fields1[hit_list_field];
                } else if (*fields1[status_field] == 'U' &&
                           *fields2[status_field] == 'C') {
                        status = "C";
                        taxid = fields2[taxid_field];
                        hit_list = fields2[hit_list_field];
                } else {
                        status = "U";
                        taxid = "0";
                        hit_list = fields1[hit_list_field];
                        total_unclassified += 1;
                }

                if (!merged_hit_list.empty()) {
                        fprintf(out, "%s\t%s\t%s\t%s\t", status,
                                fields1[header_field], taxid, fields1[len_field]);
                        write_hit_list(merged_hit_list, "", out);
                        fputc('\n', out);
                } else {
                        fprintf(out, "%s\t%s\t%s\t%s\t%s", status,
                                fields1[header_field], taxid, fields1[len_field], hit_list);
                }

                if (counters) {
                        taxid_t t = nixmans_atou64_shift(taxid, strlen(taxid));
                        taxid_t internal_taxid = taxonomy.GetInternalID(t);
                        (*counters)[internal_taxid] += kraken2::READCOUNTER(1, 0);
                }

                if (classified_headers && status[0] == 'C') {
                        fputc('>', classified_headers);
                        fwrite(fields1[header_field], strlen(fields1[header_field]),
                               1, classified_headers);
                        fputc('\n', classified_headers);
                }

                hit_list1.clear();
                hit_list2.clear();
                hit_counts.clear();
                merged_hit_list.clear();
                total_sequences += 1;
        }

        if (line1) {
                free(line1);
        }

        if (line2) {
                free(line2);
        }

        return std::tuple<size_t, size_t>(total_sequences, total_unclassified);
}

void usage(const char *prog, int err) {
        const char *errs[] = {
                "",
                "Required argument, -o, not specified",
                "Required argument, -t, not specified",
                "No input files provided",
                "Only one input file provided",
                "Too many input files provided",
                "Invalid argument",
        };

        fprintf(stderr, "Error: %s\n", errs[err]);
        fprintf(stderr, "%s -o <output> -t <merged_tax.k2d> [ -r report ] "
                "[ -c classified_headers ] input1 input2\n", prog);

        exit(err);
}

FILE *xfopen(const char *filename, const char *options) {
        FILE *f;

        if ((f = fopen(filename, options)) == nullptr) {
                err(1, "Unable to open file");
        }

        return f;
}

int main(int argc, char **argv) {

        int nflag = 0, ch;
        bool report_zeros = false;
        bool mpa_style = false;
        float confidence_threshold;
        char *report_filename = nullptr;
        char *classified_headers_filename = nullptr;
        char *merged_taxon_filename = nullptr;
        char *merged_output_filename = nullptr;
        char *input1 = nullptr;
        char *input2 = nullptr;

        nflag = 0;

        while ((ch = getopt(argc, argv, "hi:r:m:c:o:t:nz")) != -1) {
                switch (ch) {
                case 'c':
                        classified_headers_filename = optarg;
                        break;
                case 'i':
                        confidence_threshold = strtof(optarg, (char **)nullptr);
                        break;
                case 'n':
                        nflag = 1;
                        break;
                case 'o':
                        merged_output_filename = optarg;
                        break;
                case 'r':
                        report_filename = optarg;
                        break;
                case 'm':
                        report_filename = optarg;
                        mpa_style = true;
                        break;
                case 't':
                        merged_taxon_filename = optarg;
                        break;
                case 'z':
                        report_zeros = true;
                        break;
                case 'h':
                        usage(argv[0], no_error);
                default:
                        usage(argv[0], invalid_flag);
                }
        }

        argc -= optind;
        argv += optind;

        if (argc != 2) {
                switch (argc) {
                case 0:
                        usage(argv[0], no_input);
                        break;
                case 1:
                        usage(argv[0], one_input);
                        break;
                default:
                        usage(argv[0], too_much_input);
                }
        }

        if (merged_output_filename == nullptr) {
                usage(argv[0], no_output);
        }

        if (merged_taxon_filename == nullptr) {
                usage(argv[0], no_taxon);
        }

        input1 = *argv++;
        input2 = *argv++;

        FILE *l = xfopen(input1, "r");
        FILE *r = xfopen(input2, "r");
        FILE *m = xfopen(merged_output_filename, "w");
        FILE *classified_headers = nullptr;

        kraken2::taxon_counters_t *counters = nullptr;
        kraken2::Taxonomy taxonomy(merged_taxon_filename);
        taxonomy.GenerateExternalToInternalIDMap();

        if (classified_headers_filename) {
                classified_headers = xfopen(classified_headers_filename, "w");
        }
        if (report_filename) {
                counters = new kraken2::taxon_counters_t();
        }

        size_t total_seqs = 0;
        size_t total_classified = 0;
        std::tie(total_seqs, total_classified) = merge_classification_output(taxonomy, l, r, m, confidence_threshold, counters, classified_headers);

        if (report_filename != nullptr) {
                if (mpa_style) {
                        kraken2::ReportMpaStyle(report_filename, report_zeros, taxonomy, *counters);
                } else {
                        kraken2::ReportKrakenStyle(report_filename, report_zeros, false, taxonomy,
                                                   *counters, total_seqs, total_classified);
                }
        }

        delete counters;

        fclose(l);
        fclose(r);
        fclose(m);

        if (classified_headers) {
                fprintf(classified_headers, "%zu\n", total_seqs);
                fclose(classified_headers);
        }

        return 0;
}
