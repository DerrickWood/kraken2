#!/bin/bash

# Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Build a 16S database from Greengenes data

set -u  # Protect against uninitialized vars.
set -e  # Stop on error
set -o pipefail  # Stop on failures in non-final pipeline commands

FTP_SERVER="ftp://greengenes.microbio.me/"
GG_VERSION="gg_13_5"
REMOTE_DIR="$FTP_SERVER/greengenes_release/$GG_VERSION"

mkdir -p "$KRAKEN2_DB_NAME"
pushd "$KRAKEN2_DB_NAME"
mkdir -p data taxonomy library
pushd data
wget "$REMOTE_DIR/${GG_VERSION}.fasta.gz"
gunzip "${GG_VERSION}.fasta.gz"
wget "$REMOTE_DIR/${GG_VERSION}_taxonomy.txt.gz"
gunzip "${GG_VERSION}_taxonomy.txt.gz"

build_gg_taxonomy.pl "${GG_VERSION}_taxonomy.txt"
popd
mv data/names.dmp data/nodes.dmp taxonomy/
mv data/seqid2taxid.map .
mv "data/${GG_VERSION}.fasta" library/gg.fna
popd

kraken2-build --db $KRAKEN2_DB_NAME --build --threads $KRAKEN2_THREAD_CT
