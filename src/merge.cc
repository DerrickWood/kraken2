#include "kraken2_data.h"
#include "reports.h"
#include "taxonomy.h"

#include <algorithm>
#include <assert.h>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <err.h>
#include <functional>
#include <limits>
#include <sstream>
#include <string.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>
#include <stdint.h>


FILE *xfopen(const char *filename, const char *options);

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

#define BUFFER_SIZE (8 * 1024 * 1024)

class FileBuffer {
public:
        FileBuffer(const char *filename, size_t buffer_size = BUFFER_SIZE) {
                f = xfopen(filename, "r");
                if ((buffer = (char *)malloc(buffer_size)) == nullptr) {
                        errx(1, "Unable to allocate space for read buffer");
                }
                eof = 0;
                buffer_length = 0;
                buffer_capacity = buffer_size;
                buffer_position = 0;
        }

        int refill() {
                if (buffer_position > buffer_capacity) {
                        errx(1, "buffer_position (%zu) exceeds buffer_capacity (%zu)",
                             buffer_position, buffer_capacity);
                }

                size_t nread = 0;
                // we may have unread data in the buffer
                // adjust the buffer pointer and length
                // accordingly
                char *b = buffer + buffer_position;
                size_t length = buffer_capacity - buffer_position;

                while (true) {
                        size_t bytes_read = fread(b + nread, 1, length - nread, f);
                        nread += bytes_read;
                        if (nread == length || bytes_read == 0) {
                                break;
                        }
                }

                buffer_length = nread + buffer_position;
                // buffer_position = 0;

                return feof(f);
        }

        size_t copyfile(FILE *out) {
                size_t written = 0;

                do {
                        buffer_position = 0;
                        refill();
                        size_t bytes_written = 0;
                        while (true) {
                                bytes_written += fwrite(buffer + bytes_written, 1, buffer_length - bytes_written, out);
                                if (bytes_written == buffer_length) {
                                        break;
                                }
                        }
                        written += bytes_written;
                } while (buffer_length != 0);

                return written;
        }


        size_t readlines(std::vector<std::string> &lines) {
                size_t i = 0;
                size_t line_length = 0;

                while (i < lines.size()) {
                        if (buffer_position == buffer_length) {
                                buffer_position = 0;
                                eof = refill();
                        }
                        if (buffer_length == 0) {
                                i += (line_length > 0);
                                break;
                        }
                        char *newline_ptr =
                                (char *)memchr(buffer + buffer_position, '\n',
                                               buffer_length - buffer_position);
                        if (!newline_ptr) {
                                if (!eof) {
                                        // move the residual data to the beginning of the buffer
                                        memmove(
                                                buffer, &buffer[buffer_position],
                                                buffer_length - buffer_position);
                                        buffer_position =
                                            buffer_length - buffer_position;
                                        // append new data to the buffer and
                                        // reset the buffer position.
                                        eof = refill();
                                        newline_ptr = (char *)memchr(
                                                buffer + buffer_position, '\n',
                                                buffer_length - buffer_position);
                                        buffer_position = 0;
                                        if (!newline_ptr && eof) {
                                                newline_ptr = &buffer[buffer_length - 1];
                                        } else if (!newline_ptr) {
                                                // we still have not found a newline,
                                                // append the current buffer's contents
                                                // into the string and try again in
                                                // the next iteration.
                                                lines[i].resize(line_length + buffer_length);
                                                memcpy(&lines[i][line_length],
                                                       buffer, buffer_length);
                                                buffer_position = buffer_length;
                                                line_length += buffer_length;
                                                continue;
                                        }
                                } else {
                                        newline_ptr = &buffer[buffer_length - 1];
                                }
                        }

                        size_t new_length = (newline_ptr - &buffer[buffer_position]) + 1;
                        lines[i].resize(line_length + new_length);
                        memcpy(&lines[i][line_length], &buffer[buffer_position],
                               new_length);
                        buffer_position += new_length;
                        line_length = 0;
                        i += 1;
                }
                return i;
        }

        ~FileBuffer() {
                fclose(f);
        }

private:
        FILE *f;
        int eof;
        char *buffer;
        size_t buffer_length;
        size_t buffer_capacity;
        size_t buffer_position;
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

        char itoa_buffer[17] = {0};
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
                len2 = (string == nullptr ? end_p : string) - str_count - 1;

                if (!str_taxid || !str_count) {
                        break;
                }
                // count = strtol(str_count, nullptr, 10);
                // taxid = strtol(str_taxid, nullptr, 10);
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

taxid_t get_lca(kraken2::Taxonomy &taxonomy,
                std::unordered_map<uint64_t, taxid_t> &lca_cache,
                taxid_t taxid1, taxid_t taxid2) {
        if (taxid1 > taxid2) {
                taxid_t temp = taxid1;
                taxid1 = taxid2;
                taxid2 = temp;
        }
        uint64_t key = ((uint64_t)taxid1 << 32) | taxid2;
        if (lca_cache.find(key) != lca_cache.end()) {
                return lca_cache[key];
        }
        taxid1 = taxonomy.GetInternalID(taxid1);
        taxid2 = taxonomy.GetInternalID(taxid2);

        uint64_t lca = taxonomy.LowestCommonAncestor(taxid1, taxid2);

        taxid_t value = (taxid_t)taxonomy.nodes()[lca].external_id;
        lca_cache[key] = value;

        return value;
}

size_t merge_hit_lists(kraken2::Taxonomy &taxonomy,
                       std::unordered_map<uint64_t, taxid_t> &lca_cache,
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
                        final_taxid = get_lca(taxonomy, lca_cache, tc1.taxid, tc2.taxid);
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

        for (fp = fields; (*fp = strsep(&line, delim)) != nullptr;) {
                if (**fp != '\0' && ++fp >= &fields[fields_len]) {
                        break;
                }
        }
}

std::tuple<size_t, size_t> merge_classification_output_parallel(
    const char *taxonomy_filename, const char *ifn1, const char *ifn2,
    const char *ofn, const char *cfn, bool use_names,
    kraken2::taxon_counters_t *counters, float confidence_threshold,
    size_t batch_size) {

        enum {
                status_field = 0,
                header_field,
                taxid_field,
                len_field,
                hit_list_field,
        };

        FileBuffer *fr1 = new FileBuffer(ifn1);
        FileBuffer *fr2 = new FileBuffer(ifn2);
        // std::heap
        std::vector<size_t> blocks_written{std::numeric_limits<size_t>::max()};
        std::atomic<size_t> total_sequences{0};
        std::atomic<size_t> total_unclassified{0};
        std::atomic<size_t> block_count{0};
        std::atomic<bool> finished{false};

#pragma omp parallel
        {
                size_t block = 0;
                kraken2::taxon_counters_t local_counters;
                kraken2::Taxonomy taxonomy(taxonomy_filename);
                taxonomy.GenerateExternalToInternalIDMap();
                std::string out_filename = ofn;
                std::vector<std::string> lines1(batch_size);
                std::vector<std::string> lines2(batch_size);
                std::unordered_map<uint64_t, taxid_t> lca_cache;
                char *fields1[5];
                char *fields2[5];
                char scientific_name[500];

                const char *status;
                const char *taxid;

                std::vector<taxid_and_count> hit_list1;
                std::vector<taxid_and_count> hit_list2;
                std::vector<taxid_and_count> merged_hit_list;
                kraken2::taxon_counts_t hit_counts;
                char itoa_buf[17];

                size_t nlines1 = 0;
                size_t nlines2 = 0;
                ostringstream local_merge_filename;
                ostringstream local_classified_filename;

                #pragma omp single nowait
                {
                        size_t next_block = 0;
                        bool block_ready = false;
                        FILE *outfile = xfopen(ofn, "w");
                        FILE *header_outfile = nullptr;
                        if (cfn) {
                                header_outfile = xfopen(cfn, "w");
                        }
                        while (true) {
                                #pragma omp critical(heap_update)
                                {
                                        if (blocks_written.front() == next_block + 1) {
                                                std::pop_heap(blocks_written.begin(), blocks_written.end(),
                                                              std::greater<size_t>());
                                                blocks_written.pop_back();
                                                block_ready = true;
                                        }
                                }

                                if (block_ready) {
                                        next_block += 1;
                                        local_merge_filename << out_filename
                                                           << "_" << next_block;
                                                FileBuffer fb(local_merge_filename.str().c_str(),
                                                      8 * 1024);
                                        fb.copyfile(outfile);
                                        std::remove(local_merge_filename.str().c_str());
                                        local_merge_filename.str("");

                                        if (header_outfile) {
                                                local_classified_filename << cfn << "_" << next_block;
                                                FileBuffer fb2(
                                                    local_classified_filename.str().c_str(),
                                                    8 * 1024);
                                                fb2.copyfile(header_outfile);
                                                std::remove(local_classified_filename.str().c_str());
                                                local_classified_filename.str("");
                                        }

                                        block_ready = false;
                                }
                                if (block_count == next_block && finished) {
                                        break;
                                }

                        }

                        fclose(outfile);
                        if (header_outfile) {
                                fprintf(header_outfile, "%zu\n", total_sequences.load(std::memory_order_acquire));
                                fclose(header_outfile);
                        }

                }

                while (true) {
                        #pragma omp critical
                        {
                                block_count.fetch_add(1, std::memory_order_relaxed);
                                nlines1 = fr1->readlines(lines1);
                                nlines2 = fr2->readlines(lines2);
                                block = block_count;
                        }

                        local_merge_filename << out_filename << "_" << block;

                        assert(nlines1 == nlines2);

                        if (nlines1 == 0) {
                                // reduce the block count if there is no more
                                // data to be read.
                                block_count.fetch_sub(1, std::memory_order_relaxed);
                                if (counters) {
                                        #pragma omp critical(update_counters)
                                        {
                                                for (const auto &pair : local_counters) {
                                                        if (counters->find(pair.first) == counters->end()) {
                                                                (*counters)[pair.first] = pair.second;
                                                        } else {
                                                                (*counters)[pair.first] += pair.second;
                                                        }
                                                }
                                        }
                                }
                                finished = true;
                                break;
                        }

                        FILE *out = xfopen(local_merge_filename.str().c_str(), "w");
                        FILE *lcfn = nullptr;
                        if (cfn) {
                                local_classified_filename << cfn << "_" << block;
                                lcfn = xfopen(local_classified_filename.str().c_str(), "w");
                        }

                        for (size_t i = 0; i < nlines1; i++) {
                                char *line1 = &lines1[i][0];
                                char *line2 = &lines2[i][0];
                                get_fields(line1, "\t", fields1, 5);
                                get_fields(line2, "\t", fields2, 5);

                                parse_hit_list(fields1[hit_list_field], strlen(fields1[hit_list_field]), hit_list1);
                                parse_hit_list(fields2[hit_list_field], strlen(fields2[hit_list_field]), hit_list2);

                                size_t total_minimizers =
                                        merge_hit_lists(taxonomy, lca_cache, hit_counts, hit_list1, hit_list2, merged_hit_list);
                                int res = resolve_tree(taxonomy, hit_counts, total_minimizers, confidence_threshold);
                                taxid = int_to_string(res, itoa_buf);


                                if (fields1[status_field][0] == 'C' || fields2[status_field][0] == 'C') {
                                        status = "C";
                                } else {
                                        status = "U";
                                        taxid = "0";
                                        total_unclassified += 1;
                                }

                                if (use_names) {
                                        taxid_t t = nixmans_atou64_shift(taxid, strlen(taxid));
                                        taxid_t internal_taxid = taxonomy.GetInternalID(t);

                                        const char *name =
                                                taxonomy.name_data() +
                                                taxonomy.nodes()[internal_taxid].name_offset;
                                        if (!name) {
                                                name = "unclassified";
                                        }
                                        sprintf(scientific_name, "%s (taxid %s)", name, taxid);
                                        taxid = scientific_name;
                                }
                                if (!merged_hit_list.empty()) {
                                        fprintf(out, "%s\t%s\t%s\t%s\t", status,
                                                fields1[header_field], taxid, fields1[len_field]);
                                        write_hit_list(merged_hit_list, "", out);
                                        fputc('\n', out);
                                }

                                if (counters && status[0] == 'C') {
                                        taxid_t t = nixmans_atou64_shift(taxid, strlen(taxid));
                                        taxid_t internal_taxid = taxonomy.GetInternalID(t);
                                        local_counters[internal_taxid] += kraken2::READCOUNTER(1, 0);
                                }

                                if (lcfn && status[0] == 'C') {
                                        fputc('>', lcfn);
                                        fwrite(fields1[header_field], strlen(fields1[header_field]),
                                               1, lcfn);
                                        fputc('\n', lcfn);
                                }

                                hit_list1.clear();
                                hit_list2.clear();
                                hit_counts.clear();
                                merged_hit_list.clear();
                                total_sequences += 1;
                        }
                        fclose(out);
                        if (lcfn) {
                                fclose(lcfn);
                        }
#pragma omp critical(heap_update)
                        {
                                blocks_written.push_back(block);
                                std::push_heap(blocks_written.begin(), blocks_written.end(), std::greater<size_t>());
                        }
                        local_merge_filename.str("");
                        local_classified_filename.str("");
                }
        }

        delete fr1;
        delete fr2;

        return std::make_tuple<size_t, size_t>(total_sequences, total_unclassified);
}

std::tuple<size_t, size_t> merge_classification_output(
        const char *taxonomy_filename, const char *ifn1, const char *ifn2,
        const char *ofn, float confidence_threshold,
        kraken2::taxon_counters_t *counters,
        const char *classified_headers_filename, bool use_names) {

        kraken2::Taxonomy taxonomy(taxonomy_filename);
        FILE *in1 = xfopen(ifn1, "r");
        FILE *in2 = xfopen(ifn2, "r");
        FILE *out = xfopen(ofn, "w");
        FILE *classified_headers = nullptr;

        if (classified_headers_filename) {
                classified_headers = xfopen(classified_headers_filename, "w");
        }

        char *line1 = nullptr;
        char *line2 = nullptr;

        size_t line1_cap = 0;
        size_t line2_cap = 0;

        ssize_t line1_len = 0;
        ssize_t line2_len = 0;

        char *fields1[5];
        char *fields2[5];
        char scientific_name[100];

        const char *status;
        const char *taxid;

        std::unordered_map<uint64_t, taxid_t> lca_cache;

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

                parse_hit_list(fields1[hit_list_field], strlen(fields1[hit_list_field]), hit_list1);
                parse_hit_list(fields2[hit_list_field], strlen(fields2[hit_list_field]), hit_list2);

                size_t total_minimizers =
                        merge_hit_lists(taxonomy, lca_cache, hit_counts, hit_list1, hit_list2, merged_hit_list);
                int res = resolve_tree(taxonomy, hit_counts, total_minimizers, confidence_threshold);
                taxid = int_to_string(res, itoa_buf);


                if (fields1[status_field][0] == 'C' || fields2[status_field][0] == 'C') {
                        status = "C";
                } else {
                        status = "U";
                        taxid = "0";
                        total_unclassified += 1;
                }

                if (use_names) {
                        taxid_t t = nixmans_atou64_shift(taxid, strlen(taxid));
                        taxid_t internal_taxid = taxonomy.GetInternalID(t);

                        const char *name =
                            taxonomy.name_data() +
                            taxonomy.nodes()[internal_taxid].name_offset;
                        if (!name) {
                                name = "unclassified";
                        }
                        sprintf(scientific_name, "%s (taxid %s)", name, taxid);
                        taxid = scientific_name;
                }
                if (!merged_hit_list.empty()) {
                        fprintf(out, "%s\t%s\t%s\t%s\t", status,
                                fields1[header_field], taxid, fields1[len_field]);
                        write_hit_list(merged_hit_list, "", out);
                        fputc('\n', out);
                }

                if (counters && status[0] == 'C') {
                        taxid_t t = nixmans_atou64_shift(taxid, strlen(taxid));
                        taxid_t internal_taxid = taxonomy.GetInternalID(t);
                        (*counters)[internal_taxid] += kraken2::READCOUNTER(1, 0);
                }

                if (classified_headers_filename && status[0] == 'C') {
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

        fclose(in1);
        fclose(in2);
        fclose(out);

        if (classified_headers_filename) {
                fclose(classified_headers);
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
                err(1, "Unable to open file %s", filename);
        }

        return f;
}

int main(int argc, char **argv) {

        int ch;
        int threads = 1;
        int batch_size = 100;
        bool use_names = false;
        bool report_zeros = false;
        bool mpa_style = false;
        float confidence_threshold = 0.0;
        char *report_filename = nullptr;
        char *classified_headers_filename = nullptr;
        char *merged_taxon_filename = nullptr;
        char *merged_output_filename = nullptr;
        char *input1 = nullptr;
        char *input2 = nullptr;

        while ((ch = getopt(argc, argv, "b:hi:r:mc:o:t:T:nz")) != -1) {
                switch (ch) {
                case 'b':
                        batch_size = strtod(optarg, (char **)nullptr);
                        break;
                case 'c':
                        classified_headers_filename = optarg;
                        break;
                case 'i':
                        confidence_threshold = strtof(optarg, (char **)nullptr);
                        break;
                case 'n':
                        use_names = true;
                        break;
                case 'o':
                        merged_output_filename = optarg;
                        break;
                case 'r':
                        report_filename = optarg;
                        break;
                case 'm':
                        mpa_style = true;
                        break;
                case 't':
                        merged_taxon_filename = optarg;
                        break;
                case 'T':
                        threads = strtod(optarg, (char **)nullptr);
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

        kraken2::taxon_counters_t *counters = nullptr;

        if (report_filename) {
                counters = new kraken2::taxon_counters_t();
        }

        size_t total_seqs = 0;
        size_t total_unclassified = 0;
        omp_set_num_threads(threads);

        if (threads == 1) {
                std::tie(total_seqs, total_unclassified) =
                        merge_classification_output(
                                merged_taxon_filename, input1, input2, merged_output_filename,
                                confidence_threshold, counters, classified_headers_filename,
                                use_names);
        } else {

                std::tie(total_seqs, total_unclassified) =
                        merge_classification_output_parallel(
                                merged_taxon_filename, input1, input2, merged_output_filename,
                                classified_headers_filename, use_names, counters,
                                confidence_threshold, batch_size);
        }

        if (report_filename != nullptr) {
                kraken2::Taxonomy taxonomy(merged_taxon_filename);
                taxonomy.GenerateExternalToInternalIDMap();

                if (mpa_style) {
                        kraken2::ReportMpaStyle(report_filename, report_zeros, taxonomy, *counters);
                } else {
                        kraken2::ReportKrakenStyle(report_filename, report_zeros, false, taxonomy,
                                                   *counters, total_seqs, total_unclassified);
                }
        }

        // printf("total sequences: %lu\n", total_seqs);
        // printf("total_unclassified: %lu\n", total_unclassified);

        return 0;
}
