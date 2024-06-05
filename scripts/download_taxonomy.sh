#!/bin/bash

# Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Download NCBI taxonomy information for Kraken 2.
# Designed to be called by kraken2-build

set -u  # Protect against uninitialized vars.
set -e  # Stop on error

TAXONOMY_DIR="$KRAKEN2_DB_NAME/taxonomy"
NCBI_SERVER="ftp.ncbi.nlm.nih.gov"
RSYNC_SERVER="rsync://$NCBI_SERVER"
FTP_SERVER="ftp://$NCBI_SERVER"

mkdir -p "$TAXONOMY_DIR"
cd "$TAXONOMY_DIR"

function download_file() {
  file="$1"
  if [ -n "$KRAKEN2_USE_FTP" ]
  then
    wget -q ${FTP_SERVER}${file}
  else
    rsync --no-motd ${RSYNC_SERVER}${file} .
  fi
}

if [ ! -e "accmap.dlflag" ] && [ -z "$KRAKEN2_SKIP_MAPS" ]
then
  if [ -z "$KRAKEN2_PROTEIN_DB" ]
  then
    for subsection in gb wgs
    do
      1>&2 echo -n "Downloading nucleotide ${subsection} accession to taxon map..."
      download_file "/pub/taxonomy/accession2taxid/nucl_${subsection}.accession2taxid.gz"
      1>&2 echo " done."
    done
  else
    1>&2 echo -n "Downloading protein accession to taxon map..."
    download_file "/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"
    1>&2 echo " done."
  fi
  touch accmap.dlflag
  1>&2 echo "Downloaded accession to taxon map(s)"
fi

if [ ! -e "taxdump.dlflag" ]
then
  1>&2 echo -n "Downloading taxonomy tree data..."
  download_file "/pub/taxonomy/taxdump.tar.gz"
  touch taxdump.dlflag
  1>&2 echo " done."
fi

if ls | grep -q 'accession2taxid\.gz$'
then
  1>&2 echo -n "Uncompressing taxonomy data..."
  gunzip *accession2taxid.gz
  1>&2 echo " done."
fi

if [ ! -e "taxdump.untarflag" ]
then
  1>&2 echo -n "Untarring taxonomy tree data..."
  tar zxf taxdump.tar.gz
  touch taxdump.untarflag
  1>&2 echo " done."
fi
