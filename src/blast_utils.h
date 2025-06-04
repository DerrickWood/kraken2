#ifndef __BLAST_UTILS_H__
#define __BLAST_UTILS_H__

#include <inttypes.h>
#include <stdio.h>

FILE *open_file(const char *filename, const char *mode);
void *alloc_memory(void *data, uint32_t element_size, uint32_t old_size, uint32_t new_size);

uint32_t nlz(uint32_t x);
uint32_t next_power_of_2(uint32_t n);
uint32_t read_into_buffer(FILE *f, void *buffer, uint32_t element_size,
                          uint32_t buffer_len);
#endif
