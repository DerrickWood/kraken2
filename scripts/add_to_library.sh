#!/bin/bash

# Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Copy specified file into a Kraken library

set -u  # Protect against uninitialized vars.
set -e  # Stop on error

LIBRARY_DIR="$KRAKEN2_DB_NAME/library"

if [ ! -e "$1" ]
then
  1>&2 echo "Can't add \"$1\": file does not exist"
  exit 1
fi
if [ ! -f "$1" ]
then
  1>&2 echo "Can't add \"$1\": not a regular file"
  exit 1
fi

add_dir="$LIBRARY_DIR/added"
mkdir -p "$add_dir"
prelim_map=$(cp_into_tempfile.pl -t "prelim_map_XXXXXXXXXX" -d "$add_dir" -s txt /dev/null)
scan_fasta_file.pl "$1" > "$prelim_map"

filename=$(cp_into_tempfile.pl -t "XXXXXXXXXX" -d "$add_dir" -s fna "$1")

if [ -n "$KRAKEN2_MASK_LC" ]; then
  1>&2 echo -n "Masking low-complexity regions of new file..."
  mask_low_complexity.sh $filename
  1>&2 echo " done."
fi

1>&2 echo "Added \"$1\" to library ($KRAKEN2_DB_NAME)"
