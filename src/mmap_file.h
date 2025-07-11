/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_MMAP_FILE_H_
#define KRAKEN2_MMAP_FILE_H_

#include "kraken2_headers.h"

namespace kraken2 {
  class MMapFile {
    public:
      MMapFile();
      MMapFile(MMapFile &&o) {
        valid_ = o.valid_;
        o.valid_ = false;
        fd_ = o.fd_;
        o.fd_ = -1;
        fptr_ = o.fptr_;
        o.fptr_ = NULL;
        filesize_ = o.filesize_;
        o.filesize_ = 0;
      };
    ~MMapFile();
    void OpenFile(const std::string &filename, int mode = O_RDONLY, int map_flags = -1, int prot_flags = -1, size_t size = 0);
    void OpenFile(const char *filename, int mode = O_RDONLY, int map_flags = -1, int prot_flags = -1, size_t size = 0);

    char *fptr();
    size_t filesize();

    void LoadFile();
    void SyncFile();
    void CloseFile();

    private:
    MMapFile(const MMapFile &rhs);
    MMapFile& operator=(const MMapFile &rhs);

    bool valid_;
    int fd_;
    char *fptr_;
    size_t filesize_;
  };
}

#endif
