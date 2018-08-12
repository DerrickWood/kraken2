#!/bin/bash

# Copyright 2013-2018, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Masks low complexity sequences in the database using the dustmasker (nucl.)
# or segmasker (prot.) programs from NCBI.

set -u  # Protect against uninitialized vars.
set -e  # Stop on error
set -o pipefail  # Stop on failures in non-final pipeline commands

target="$1"

MASKER="dustmasker"
if [ -n "$KRAKEN2_PROTEIN_DB" ]; then
  MASKER="segmasker"
fi

if ! which $MASKER > /dev/null; then
  echo "Unable to find $MASKER in path, can't mask low-complexity sequences"
  exit 1
fi

function mask_data_chunk () {
  # Removes empty records and performs masking, all in pipes
  MASKER=$1
  awk -v RS=">" -v FS="\n" -v ORS="" ' { if ($2) print ">"$0 } ' |\
  $MASKER -in - -outfmt fasta |\
  sed -e '/^>/!s/[a-z]/x/g'
}
export -f mask_data_chunk

if [ -d $target ]; then
  for file in $(find $target '(' -name '*.fna' -o -name '*.faa' ')'); do
    if [ ! -e "$file.masked" ]; then
      cat $file | parallel --pipe --recstart '>' --blocksize 100M mask_data_chunk $MASKER > "$file.tmp"
      mv "$file.tmp" $file
      touch "$file.masked"
    fi
  done
elif [ -f $target ]; then
  if [ ! -e "$target.masked" ]; then
    cat $target | parallel --pipe --recstart '>' --blocksize 100M mask_data_chunk $MASKER > "$target.tmp"
    mv "$target.tmp" $target
    touch "$target.masked"
  fi
else
  echo "Target $target must be directory or regular file, aborting masking"
  exit 1
fi
