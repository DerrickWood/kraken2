#!/bin/bash

# Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Masks low complexity sequences in the database using the dustmasker (nucl.)
# or segmasker (prot.) programs from NCBI.

set -u  # Protect against uninitialized vars.
set -e  # Stop on error
set -o pipefail  # Stop on failures in non-final pipeline commands

target="$1"

MASKER="k2mask"
if [ -n "$KRAKEN2_PROTEIN_DB" ]; then
  MASKER="segmasker"
fi

if ! which $MASKER > /dev/null; then
  echo "Unable to find $MASKER in path, can't mask low-complexity sequences"
  exit 1
fi

if [ -d $target ]; then
  for file in $(find $target '(' -name '*.fna' -o -name '*.faa' ')'); do
    if [ ! -e "$file.masked" ]; then
      $MASKER -in $file -outfmt fasta | sed -e '/^>/!s/[a-z]/x/g' > "$file.tmp"
      mv "$file.tmp" $file
      touch "$file.masked"
    fi
  done
elif [ -f $target ]; then
  if [ ! -e "$target.masked" ]; then
    $MASKER -in $target -outfmt fasta | sed -e '/^>/!s/[a-z]/x/g' > "$target.tmp"
    mv "$target.tmp" $target
    touch "$target.masked"
  fi
else
  echo "Target $target must be directory or regular file, aborting masking"
  exit 1
fi
