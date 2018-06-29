#!/bin/bash

# Copyright 2013-2018, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Build a 16S database from Silva data

set -u  # Protect against uninitialized vars.
set -e  # Stop on error
set -o pipefail  # Stop on failures in non-final pipeline commands

FTP_SERVER="ftp://ftp.arb-silva.de/"
SILVA_VERSION="132"
REMOTE_DIR="$FTP_SERVER/release_$SILVA_VERSION/Exports"
FASTA_FILENAME="SILVA_${SILVA_VERSION}_SSURef_Nr99_tax_silva.fasta"
TAXO_PREFIX="tax_slv_ssu_$SILVA_VERSION"

mkdir -p "$KRAKEN2_DB_NAME"
pushd "$KRAKEN2_DB_NAME"
mkdir -p data taxonomy library
pushd data
wget "$REMOTE_DIR/${FASTA_FILENAME}.gz"
gunzip "${FASTA_FILENAME}.gz"
wget "$REMOTE_DIR/taxonomy/${TAXO_PREFIX}.acc_taxid"
wget "$REMOTE_DIR/taxonomy/${TAXO_PREFIX}.txt"

build_silva_taxonomy.pl "${TAXO_PREFIX}.txt"
popd
mv data/names.dmp data/nodes.dmp taxonomy/
mv data/${TAXO_PREFIX}.acc_taxid seqid2taxid.map
sed -e '/^>/!y/U/T/' "data/$FASTA_FILENAME" > library/silva.fna
popd

kraken2-build --db $KRAKEN2_DB_NAME --build --threads $KRAKEN2_THREAD_CT
