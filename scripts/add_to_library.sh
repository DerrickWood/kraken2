#!/bin/bash

# Copyright 2013-2018, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Copy specified file into a Kraken library

set -u  # Protect against uninitialized vars.
set -e  # Stop on error

LIBRARY_DIR="$KRAKEN2_DB_NAME/library"

if [ ! -e "$1" ]
then
  echo "Can't add \"$1\": file does not exist"
  exit 1
fi
if [ ! -f "$1" ]
then
  echo "Can't add \"$1\": not a regular file"
  exit 1
fi

add_dir="$LIBRARY_DIR/added"
mkdir -p "$add_dir"
scan_fasta_file.pl "$1" > "$add_dir/temp_map.txt"

filename=$(cp_into_tempfile.pl -t "XXXXXXXXXX" -d "$add_dir" -s fna "$1")

if [ -n "$KRAKEN2_MASK_LC" ]; then
  echo -n "Masking low-complexity regions of new file..."
  mask_low_complexity.sh $filename
  echo "done."
fi

cat "$add_dir/temp_map.txt" >> "$add_dir/prelim_map.txt"
rm "$add_dir/temp_map.txt"

echo "Added \"$1\" to library ($KRAKEN2_DB_NAME)"
