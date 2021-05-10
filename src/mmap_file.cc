/*
 * Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include "mmap_file.h"

using std::string;

namespace kraken2 {

MMapFile::MMapFile() {
  valid_ = false;
  fptr_ = NULL;
  filesize_ = 0;
  fd_ = -1;
}

void MMapFile::OpenFile(const string &filename_str, int mode, int prot_flags,
  int map_flags, size_t size)
{
  const char *filename = filename_str.c_str();
  OpenFile(filename, mode, map_flags, prot_flags, size);
}

void MMapFile::OpenFile(const char *filename, int mode, int prot_flags,
  int map_flags, size_t size)
{
  if (mode & O_APPEND || (mode & O_ACCMODE) == O_WRONLY)
    errx(EX_SOFTWARE, "illegal mode passed to MMapFile");
  if (prot_flags < 0)
    prot_flags = (mode & O_ACCMODE) == O_RDONLY
                  ? PROT_READ
                  : PROT_READ | PROT_WRITE;
  if (map_flags < 0)
    map_flags = (mode & O_ACCMODE) == O_RDONLY ? MAP_PRIVATE : MAP_SHARED;

  fd_ = open(filename, mode, 0666);
  if (fd_ < 0)
    err(EX_OSERR, "unable to open %s", filename);

  if (mode & O_CREAT) {
    if (lseek(fd_, size - 1, SEEK_SET) < 0)
      err(EX_OSERR, "unable to lseek (%s)", filename);
    if (write(fd_, "", 1) < 0)
      err(EX_OSERR, "write error (%s)", filename);
    filesize_ = size;
  }
  else {
    struct stat sb;
    if (fstat(fd_, &sb) < 0)
      err(EX_OSERR, "unable to fstat %s", filename);
    filesize_ = sb.st_size;
  }

  fptr_ = (char *) mmap(0, filesize_, prot_flags, map_flags, fd_, 0);
  if (fptr_ == MAP_FAILED) {
    err(EX_OSERR, "unable to mmap %s", filename);
  }
  valid_ = true;
}

// Basically a cat operation, loads file into OS cache
// I don't use MAP_POPULATE to do this because of portability issues
void MMapFile::LoadFile() {
  int thread_ct = 1;
  int thread = 0;
  #ifdef _OPENMP
  int old_thread_ct = omp_get_max_threads();
  // Don't use more than 4 threads, don't want to hammer disk
  if (old_thread_ct > 4)
    omp_set_num_threads(4);
  thread_ct = omp_get_max_threads();
  #endif

  size_t page_size = getpagesize();
  char buf[thread_ct][page_size];

  #pragma omp parallel
  {
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif

    #pragma omp for schedule(dynamic)
    for (size_t pos = 0; pos < filesize_; pos += page_size) {
      size_t this_page_size = filesize_ - pos;
      if (this_page_size > page_size)
        this_page_size = page_size;
      memcpy(buf[thread], fptr_ + pos, this_page_size);
    }
  }

  #ifdef _OPENMP
  omp_set_num_threads(old_thread_ct);
  #endif
}

char * MMapFile::fptr() {
  return valid_ ? fptr_ : NULL;
}

size_t MMapFile::filesize() {
  return valid_ ? filesize_ : 0;
}

MMapFile::~MMapFile() {
  CloseFile();
}

void MMapFile::SyncFile() {
  msync(fptr_, filesize_, MS_SYNC);
}

void MMapFile::CloseFile() {
  if (! valid_)
    return;
  SyncFile();
  munmap(fptr_, filesize_);
  close(fd_);
  valid_ = false;
}

} // namespace
