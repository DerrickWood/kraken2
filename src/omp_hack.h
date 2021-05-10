/*
 * Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_OMP_HACK_H_
#define KRAKEN2_OMP_HACK_H_

#ifdef _OPENMP
#include <omp.h>
#else

// This file defines a datatype and declares functions supplied by the
// omp.h file to allow compilation on systems that do not have OpenMP
// supported by the installed GCC.  Note that these functions are
// effectively no-op functions that will allow single-threaded operation only.

typedef int omp_lock_t;

int omp_get_thread_num();
int omp_get_max_threads();
void omp_set_num_threads(int num);
void omp_init_lock(omp_lock_t *lock);
void omp_destroy_lock(omp_lock_t *lock);
void omp_set_lock(omp_lock_t *lock);
void omp_unset_lock(omp_lock_t *lock);
int omp_test_lock(omp_lock_t *lock);

#endif  // defined _OPENMP

#endif
