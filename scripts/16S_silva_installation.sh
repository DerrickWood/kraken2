#!/bin/bash

# Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Build a 16S database from Silva data

set -u  # Protect against uninitialized vars.
set -e  # Stop on error
set -o pipefail  # Stop on failures in non-final pipeline commands

HTTPS_SERVER="https://ftp.arb-silva.de/"
SILVA_VERSION="138_1"
REMOTE_DIR="$HTTPS_SERVER/release_$SILVA_VERSION/Exports"
FASTA_FILENAME="SILVA_${SILVA_VERSION/_/.}_SSURef_NR99_tax_silva.fasta"
TAXO_PREFIX="tax_slv_ssu_${SILVA_VERSION/_/.}"

mkdir -p "$KRAKEN2_DB_NAME"
pushd "$KRAKEN2_DB_NAME"
mkdir -p data taxonomy library
pushd data
wget "$REMOTE_DIR/${FASTA_FILENAME}.gz"
gunzip "${FASTA_FILENAME}.gz"
wget "$REMOTE_DIR/taxonomy/${TAXO_PREFIX}.acc_taxid.gz"
gunzip "${TAXO_PREFIX}.acc_taxid.gz"
wget "$REMOTE_DIR/taxonomy/${TAXO_PREFIX}.txt.gz"
gunzip "${TAXO_PREFIX}.txt.gz"

build_silva_taxonomy.pl "${TAXO_PREFIX}.txt"
popd
mv data/names.dmp data/nodes.dmp taxonomy/
mv data/${TAXO_PREFIX}.acc_taxid seqid2taxid.map
sed -e '/^>/!y/U/T/' "data/$FASTA_FILENAME" > library/silva.fna
popd

kraken2-build --db $KRAKEN2_DB_NAME --build --threads $KRAKEN2_THREAD_CT
