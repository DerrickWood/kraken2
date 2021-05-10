#!/bin/bash

# Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Build a 16S database from RDP data

set -u  # Protect against uninitialized vars.
set -e  # Stop on error
set -o pipefail  # Stop on failures in non-final pipeline commands

HTTP_SERVER="http://rdp.cme.msu.edu/"
REMOTE_DIR="$HTTP_SERVER/download/"

mkdir -p "$KRAKEN2_DB_NAME"
pushd "$KRAKEN2_DB_NAME"
mkdir -p data taxonomy library
pushd data
wget "$REMOTE_DIR/current_Bacteria_unaligned.fa.gz"
gunzip "current_Bacteria_unaligned.fa.gz"
wget "$REMOTE_DIR/current_Archaea_unaligned.fa.gz"
gunzip "current_Archaea_unaligned.fa.gz"

build_rdp_taxonomy.pl current_*_unaligned.fa
popd
mv data/names.dmp data/nodes.dmp taxonomy/
mv data/seqid2taxid.map .
for file in data/*.fa; do
  mv $file library/$(basename $file .fa).fna
done
popd

kraken2-build --db $KRAKEN2_DB_NAME --build --threads $KRAKEN2_THREAD_CT
