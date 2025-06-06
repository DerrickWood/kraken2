CXX ?= g++
KRAKEN2_SKIP_FOPENMP ?= -fopenmp
CXXFLAGS = $(KRAKEN2_SKIP_FOPENMP) -Wall -std=c++11 -O3
CXXFLAGS += -DLINEAR_PROBING
CFLAGS = -Wall -std=c99 -O2 -g3

.PHONY: all clean install

PROGS = estimate_capacity build_db classify dump_table lookup_accession_numbers k2mask blast_to_fasta

all: $(PROGS)

install: $(PROGS)
	cp $(PROGS) "$(KRAKEN2_DIR)/"

clean:
	rm -f *.o $(PROGS)

taxonomy.o: taxonomy.cc taxonomy.h mmap_file.h
hyperloglogplus.o: hyperloglogplus.cc hyperloglogplus.h
mmap_file.o: mmap_file.cc mmap_file.h
compact_hash.o: compact_hash.cc compact_hash.h kv_store.h mmap_file.h kraken2_data.h
mmscanner.o: mmscanner.cc mmscanner.h
seqreader.o: seqreader.cc seqreader.h kseq.h
omp_hack.o: omp_hack.cc omp_hack.h
reports.o: reports.cc reports.h kraken2_data.h
aa_translate.o: aa_translate.cc aa_translate.h
utilities.o: utilities.cc utilities.h

classify.o: classify.cc kraken2_data.h kv_store.h taxonomy.h seqreader.h mmscanner.h compact_hash.h aa_translate.h reports.h utilities.h readcounts.h
dump_table.o: dump_table.cc compact_hash.h taxonomy.h mmscanner.h kraken2_data.h reports.h
estimate_capacity.o: estimate_capacity.cc kv_store.h mmscanner.h seqreader.h utilities.h
build_db.o: build_db.cc taxonomy.h mmscanner.h seqreader.h compact_hash.h kv_store.h kraken2_data.h utilities.h
lookup_accession_numbers.o: lookup_accession_numbers.cc mmap_file.h utilities.h
k2mask.o: k2mask.cc gzstream.h threadpool.h

blast_defline.o: blast_defline.c
blast_to_fasta.o: blast_to_fasta.c
blast_utils.o: blast_utils.c

build_db: build_db.o mmap_file.o compact_hash.o taxonomy.o seqreader.o mmscanner.o omp_hack.o utilities.o
	$(CXX) $(CXXFLAGS) -o $@ build_db.o mmap_file.o compact_hash.o taxonomy.o seqreader.o mmscanner.o omp_hack.o utilities.o

classify: classify.o reports.o hyperloglogplus.o mmap_file.o compact_hash.o taxonomy.o seqreader.o mmscanner.o omp_hack.o aa_translate.o utilities.o
	$(CXX) $(CXXFLAGS) -o $@ classify.o reports.o hyperloglogplus.o mmap_file.o compact_hash.o taxonomy.o seqreader.o mmscanner.o omp_hack.o aa_translate.o utilities.o

estimate_capacity: estimate_capacity.o seqreader.o mmscanner.o omp_hack.o utilities.o
	$(CXX) $(CXXFLAGS) -o $@ estimate_capacity.o seqreader.o mmscanner.o omp_hack.o utilities.o

dump_table: dump_table.o mmap_file.o compact_hash.o omp_hack.o taxonomy.o reports.o hyperloglogplus.o
	$(CXX) $(CXXFLAGS) -o $@ dump_table.o mmap_file.o compact_hash.o omp_hack.o taxonomy.o reports.o hyperloglogplus.o

lookup_accession_numbers: lookup_accession_numbers.o mmap_file.o omp_hack.o utilities.o
	$(CXX) $(CXXFLAGS) -o $@ lookup_accession_numbers.o mmap_file.o omp_hack.o utilities.o

k2mask: k2mask.o seqreader.o
	$(CXX) $(CXXFLAGS) -o $@ k2mask.o seqreader.o -lz -pthread

blast_to_fasta: blast_to_fasta.o blast_defline.o blast_utils.o
	$(CC) $(CFLAGS) -o $@  blast_to_fasta.o blast_utils.o blast_defline.o -lm -lc #-fsanitize=address,undefined
