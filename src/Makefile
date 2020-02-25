CXXFLAGS = -fopenmp -Wall -std=c++11 -O3
CXXFLAGS += -DLINEAR_PROBING

.PHONY: all clean install

PROGS = estimate_capacity build_db classify dump_table lookup_accession_numbers

all: $(PROGS)

install: $(PROGS)
	cp $(PROGS) $(KRAKEN2_DIR)/

clean:
	rm -f *.o $(PROGS)

taxonomy.o: taxonomy.cc taxonomy.h mmap_file.h
hyperloglogplus.o: hyperloglogplus.cc hyperloglogplus.h 
mmap_file.o: mmap_file.cc mmap_file.h
compact_hash.o: compact_hash.cc compact_hash.h kv_store.h mmap_file.h
mmscanner.o: mmscanner.cc mmscanner.h
seqreader.o: seqreader.cc seqreader.h
omp_hack.o: omp_hack.cc omp_hack.h
reports.o: reports.cc reports.h
aa_translate.o: aa_translate.cc aa_translate.h
utilities.o: utilities.cc utilities.h

build_db: build_db.cc mmap_file.o compact_hash.o taxonomy.o seqreader.o mmscanner.o omp_hack.o utilities.o

classify: classify.cc reports.o hyperloglogplus.o mmap_file.o compact_hash.o taxonomy.o seqreader.o mmscanner.o omp_hack.o aa_translate.o utilities.o

estimate_capacity: estimate_capacity.cc seqreader.o mmscanner.o omp_hack.o utilities.o

dump_table: dump_table.cc mmap_file.o compact_hash.o omp_hack.o taxonomy.o reports.o

lookup_accession_numbers: lookup_accession_numbers.cc mmap_file.o omp_hack.o utilities.o
