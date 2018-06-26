#!/bin/bash

# Copyright 2013-2018, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Build a Kraken 2 database
# Designed to be called by kraken2-build

set -u  # Protect against uninitialized vars.
set -e  # Stop on error
set -o pipefail  # Stop on failures in non-final pipeline commands

function finalize_file() {
  mv $1.tmp $1
}

function get_current_time() {
  date "+%s.%N"
}

function report_time_elapsed() {
  curr_time=$(get_current_time)
  perl -e '$time = $ARGV[1] - $ARGV[0];' \
       -e '$sec = int($time); $nsec = $time - $sec;' \
       -e '$min = int($sec/60); $sec %= 60;' \
       -e '$hr = int($min/60); $min %= 60;' \
       -e 'print "${hr}h" if $hr;' \
       -e 'print "${min}m" if $min || $hr;' \
       -e 'printf "%.3fs", $sec + $nsec;' \
       $1 $curr_time
}

function list_sequence_files() {
  find library/ '(' -name '*.fna' -o -name '*.faa' ')' -print0
}

start_time=$(get_current_time)

DATABASE_DIR="$KRAKEN2_DB_NAME"

if [ ! -d "$DATABASE_DIR" ]
then
  echo "Can't find Kraken 2 DB directory \"$KRAKEN2_DB_NAME\""
  exit 1
fi
cd "$DATABASE_DIR"

if [ ! -d "taxonomy/" ]
then
  echo "Can't find taxonomy/ subdirectory in database directory, exiting."
  exit 1
fi

if [ ! -d "library/" ]
then
  echo "Can't find library/ subdirectory in database directory, exiting."
  exit 1
fi

KRAKEN2XFLAG=""
if [ -n "$KRAKEN2_PROTEIN_DB" ]
then
  KRAKEN2XFLAG="-X"
fi

echo "Creating sequence ID to taxonomy ID map (step 1)..."
seqid2taxid_map_file=seqid2taxid.map
if [ -e "$seqid2taxid_map_file" ]; then
  echo "Sequence ID to taxonomy ID map already present, skipping map creation."
else
  step_time=$(get_current_time)
  find library/ -maxdepth 2 -name prelim_map.txt | xargs cat > taxonomy/prelim_map.txt
  if [ ! -s "taxonomy/prelim_map.txt" ]; then
    echo "No preliminary seqid/taxid mapping files found, aborting."
    exit 1
  fi
  grep "^TAXID" taxonomy/prelim_map.txt | cut -f 2- > $seqid2taxid_map_file.tmp || true
  if grep "^ACCNUM" taxonomy/prelim_map.txt | cut -f 2- > accmap_file.tmp; then
    if compgen -G "taxonomy/*.accession2taxid" > /dev/null; then
      lookup_accession_numbers.pl accmap_file.tmp taxonomy/*.accession2taxid > seqid2taxid_acc.tmp
      cat seqid2taxid_acc.tmp >> $seqid2taxid_map_file.tmp
      rm seqid2taxid_acc.tmp
    else
      echo "Accession to taxid map files are required to build this DB."
      echo "Run 'kraken2-build --db $KRAKEN2_DB_NAME --download-taxonomy' again?"
      exit 1
    fi
  fi
  rm -f accmap_file.tmp
  finalize_file $seqid2taxid_map_file
  echo "Sequence ID to taxonomy ID map complete. [$(report_time_elapsed $step_time)]"
fi

echo "Estimating required capacity (step 2)..."

step_time=$(get_current_time)
estimate=$(list_sequence_files | xargs -0 cat | estimate_capacity -k $KRAKEN2_KMER_LEN -l $KRAKEN2_MINIMIZER_LEN -S $KRAKEN2_SEED_TEMPLATE -p $KRAKEN2_THREAD_CT $KRAKEN2XFLAG )
required_capacity=$(perl -le 'print int(shift() / 0.7)' $estimate);

echo "Estimated capacity requirement: $(( required_capacity * 4 )) bytes"
echo "Capacity estimation complete. [$(report_time_elapsed $step_time)]"

echo "Building database files (step 3)..."

if [ -e "hash.k2d" ]
then
  echo "Hash table already present, skipping database file build."
else
  step_time=$(get_current_time)
  list_sequence_files | xargs -0 cat | \
    build_db -k $KRAKEN2_KMER_LEN -l $KRAKEN2_MINIMIZER_LEN -S $KRAKEN2_SEED_TEMPLATE $KRAKEN2XFLAG \
             -H hash.k2d.tmp -t taxo.k2d.tmp -o opts.k2d.tmp -n taxonomy/ -m $seqid2taxid_map_file \
             -c $required_capacity -p $KRAKEN2_THREAD_CT
  finalize_file taxo.k2d
  finalize_file opts.k2d
  finalize_file hash.k2d
  echo "Database files completed. [$(report_time_elapsed $step_time)]"
fi

echo "Database construction complete. [Total: $(report_time_elapsed $start_time)]"
