#include <err.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
// Taken from Hackers Delight
uint32_t nlz(uint32_t x) {
  uint32_t n;

  if (x == 0) {
    return 32;
  }
  n = 1;
  if ((x >> 16) == 0) {
    n = n + 16;
    x = x << 16;
  }
  if ((x >> 24) == 0) {
    n = n + 8;
    x = x << 8;
  }
  if ((x >> 28) == 0) {
    n = n + 4;
    x = x << 4;
  }
  if ((x >> 30) == 0) {
    n = n + 2;
    x = x << 2;
  }
  n = n - (x >> 31);

  return n;
}

uint32_t next_power_of_2(uint32_t n) {
  return 1 << (32 - nlz(n - 1));
}

void *alloc_memory(void *data, uint32_t element_size, uint32_t old_size,
                   uint32_t new_size) {
  void *alloc;

  if (old_size == 0 && new_size > 0) {
    alloc = malloc(element_size * new_size);
  } else if (new_size > old_size) {
    alloc = realloc(data, new_size * element_size);
  } else {
    alloc = data;
  }
  /* if (new_size > old_size) { */
  /*   bzero(alloc, new_size - old_size); */
  /* } */

  return alloc;
}

FILE *open_file(const char *filename, const char *mode) {
  if (filename == NULL) {
    return NULL;
  }

  FILE *f = fopen(filename, mode);
  if (f == NULL) {
    errx(1, "fopen");
  }

  return f;
}

uint32_t read_into_buffer(FILE *f, void *buffer, uint32_t element_size,
                          uint32_t buffer_len) {
  /* uint32_t bytes_read = 0; */
  /* uint32_t total_size = buffer_len * element_size; */
  /* uint32_t bytes_left = buffer_len; */

  /* do { */
  /*   bytes_read = fread(buffer, element_size, bytes_left, f); */
  /*   if (bytes_read == -1) { */
  /*     errx(1, "fread"); */
  /*   } */
  /*   bytes_left -= bytes_read; */
  /* } while (bytes_read == 0 || bytes_left != 0); */

  /* return total_size - bytes_left; */

  return fread(buffer, element_size, buffer_len, f);
}
