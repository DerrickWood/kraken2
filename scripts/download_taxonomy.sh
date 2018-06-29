#!/bin/bash

# Copyright 2013-2018, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Download NCBI taxonomy information for Kraken 2.
# Designed to be called by kraken2-build

set -u  # Protect against uninitialized vars.
set -e  # Stop on error

TAXONOMY_DIR="$KRAKEN2_DB_NAME/taxonomy"
NCBI_SERVER="ftp.ncbi.nlm.nih.gov"
FTP_SERVER="ftp://$NCBI_SERVER"

mkdir -p "$TAXONOMY_DIR"
cd "$TAXONOMY_DIR"

function wget_clobber() {
  filename=$(basename $1)
  rm -f $filename
  wget $1
}

function check_and_download_nucl_accmap_file() {
  flag_file="accmap_${1}.dlflag"
  if [ ! -e "$flag_file" ]
  then
    wget_clobber $FTP_SERVER/pub/taxonomy/accession2taxid/nucl_${1}.accession2taxid.gz
    touch "$flag_file"
  fi
}

if [ ! -e "accmap.dlflag" ]
then
  if [ -z "$KRAKEN2_PROTEIN_DB" ]
  then
    check_and_download_nucl_accmap_file "est"
    check_and_download_nucl_accmap_file "gb"
    check_and_download_nucl_accmap_file "gss"
    check_and_download_nucl_accmap_file "wgs"
  else
    wget_clobber $FTP_SERVER/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
  fi
  touch accmap.dlflag
  echo "Downloaded accession to taxon map(s)"
fi

if [ ! -e "taxdump.dlflag" ]
then
  wget_clobber $FTP_SERVER/pub/taxonomy/taxdump.tar.gz
  touch taxdump.dlflag
  echo "Downloaded taxonomy tree data"
fi

if ls | grep -q 'accession2taxid\.gz$'
then
  echo -n "Uncompressing taxonomy data... "
  gunzip *accession2taxid.gz
  echo "done."
fi

if [ ! -e "taxdump.untarflag" ]
then
  echo -n "Untarring taxonomy tree data... "
  tar zxf taxdump.tar.gz
  touch taxdump.untarflag
  echo "done."
fi
