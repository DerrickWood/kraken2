add_executable(build_db
        build_db.cc
        mmap_file.cc
        compact_hash.cc
        taxonomy.cc
        seqreader.cc
        mmscanner.cc
        omp_hack.cc
        utilities.cc)

add_executable(classify
        classify.cc
        reports.cc
        mmap_file.cc
        compact_hash.cc
        taxonomy.cc
        seqreader.cc
        mmscanner.cc
        omp_hack.cc
        aa_translate.cc
        utilities.cc
        hyperloglogplus.cc)

add_executable(estimate_capacity
        estimate_capacity.cc
        seqreader.cc
        mmscanner.cc
        omp_hack.cc
        utilities.cc)

add_executable(dump_table
        dump_table.cc
        mmap_file.cc
        compact_hash.cc
        omp_hack.cc
        taxonomy.cc
        reports.cc
        hyperloglogplus.cc)

add_executable(lookup_accession_numbers
        lookup_accession_numbers.cc
        mmap_file.cc
        omp_hack.cc
        utilities.cc)
