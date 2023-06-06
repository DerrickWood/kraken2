#!/bin/bash

# Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Download specific genomic libraries for use with Kraken 2.
# Supported libraries were chosen based on support from NCBI's FTP site
#   in easily obtaining a good collection of genomic data.  Others may
#   be added upon popular demand.

set -u  # Protect against uninitialized vars.
set -e  # Stop on error

LIBRARY_DIR="$KRAKEN2_DB_NAME/library"
NCBI_SERVER="ftp.ncbi.nlm.nih.gov"
FTP_SERVER="ftp://$NCBI_SERVER"
RSYNC_SERVER="rsync://$NCBI_SERVER"
THIS_DIR=$PWD

library_name="$1"
ftp_subdir=$library_name
library_file="library.fna"
if [ -n "$KRAKEN2_PROTEIN_DB" ]; then
  library_file="library.faa"
fi

function download_file() {
  file="$1"
  if [ -n "$KRAKEN2_USE_FTP" ]
  then
    wget -q ${FTP_SERVER}${file}
  else
    rsync --no-motd ${RSYNC_SERVER}${file} .
  fi
}

case $library_name in
  "archaea" | "bacteria" | "viral" | "fungi" | "plant" | "human" | "protozoa")
    mkdir -p $LIBRARY_DIR/$library_name
    cd $LIBRARY_DIR/$library_name
    rm -f assembly_summary.txt
    remote_dir_name=$library_name
    if [ "$library_name" = "human" ]; then
      remote_dir_name="vertebrate_mammalian/Homo_sapiens"
    fi
    if ! download_file "/genomes/refseq/$remote_dir_name/assembly_summary.txt"; then
      1>&2 echo "Error downloading assembly summary file for $library_name, exiting."
      exit 1
    fi
    if [ "$library_name" = "human" ]; then
      grep "Genome Reference Consortium" assembly_summary.txt > x
      mv x assembly_summary.txt
    fi
    rm -rf all/ library.f* manifest.txt rsync.err
    rsync_from_ncbi.pl assembly_summary.txt
    scan_fasta_file.pl $library_file >> prelim_map.txt
    ;;
  "plasmid")
    mkdir -p $LIBRARY_DIR/plasmid
    cd $LIBRARY_DIR/plasmid
    rm -f library.f* plasmid.*
    ## This is staying FTP only D/L for now
    1>&2 echo -n "Downloading plasmid files from FTP..."
    wget -q --no-remove-listing --spider $FTP_SERVER/genomes/refseq/plasmid/
    if [ -n "$KRAKEN2_PROTEIN_DB" ]; then
      awk '{ print $NF }' .listing | perl -ple 'tr/\r//d' | grep '\.faa\.gz' > manifest.txt
    else
      awk '{ print $NF }' .listing | perl -ple 'tr/\r//d' | grep '\.fna\.gz' > manifest.txt
    fi
    cat manifest.txt | xargs -n1 -I{} wget -q $FTP_SERVER/genomes/refseq/plasmid/{}
    cat manifest.txt | xargs -n1 -I{} gunzip -c {} > $library_file
    rm -f plasmid.* .listing
    scan_fasta_file.pl $library_file > prelim_map.txt
    1>&2 echo " done."
    ;;
  "nr" | "nt")
    protein_lib=0
    if [ "$library_name" = "nr" ]; then
      protein_lib=1
    fi
    if (( protein_lib == 1 )) && [ -z "$KRAKEN2_PROTEIN_DB" ]; then
      1>&2 echo "$library_name is a protein database, and the Kraken DB specified is nucleotide"
      exit 1
    fi
    mkdir -p $LIBRARY_DIR/$library_name
    cd $LIBRARY_DIR/$library_name
    rm -f $library_name.gz
    1>&2 echo -n "Downloading $library_name database from server... "
    download_file "/blast/db/FASTA/$library_name.gz"
    1>&2 echo "done."
    1>&2 echo -n "Uncompressing $library_name database..."
    gunzip $library_name.gz
    mv $library_name $library_file
    1>&2 echo "done."
    1>&2 echo -n "Parsing $library_name FASTA file..."
    # The nr/nt files tend to have non-standard sequence IDs, so
    # --lenient is used here.
    scan_fasta_file.pl --lenient $library_file >> prelim_map.txt
    1>&2 echo "done."
    ;;
  "UniVec" | "UniVec_Core")
    if [ -n "$KRAKEN2_PROTEIN_DB" ]; then
      1>&2 echo "$library_name is for nucleotide databases only"
      exit 1
    fi
    mkdir -p $LIBRARY_DIR/$library_name
    cd $LIBRARY_DIR/$library_name
    1>&2 echo -n "Downloading $library_name data from server... "
    download_file "/pub/UniVec/$library_name"
    1>&2 echo "done."
    # 28384: "other sequences"
    special_taxid=28384
    1>&2 echo -n "Adding taxonomy ID of $special_taxid to all sequences... "
    sed -e "s/^>/>kraken:taxid|$special_taxid|/" $library_name > library.fna
    scan_fasta_file.pl library.fna > prelim_map.txt
    1>&2 echo "done."
    ;;
  *)
    1>&2 echo "Unsupported library.  Valid options are: "
    1>&2 echo "  archaea bacteria viral fungi plant protozoa human plasmid"
    1>&2 echo "  nr nt UniVec UniVec_Core"
    exit 1
    ;;
esac

if [ -n "$KRAKEN2_MASK_LC" ]; then
  1>&2 echo -n "Masking low-complexity regions of downloaded library..."
  mask_low_complexity.sh .
  1>&2 echo " done."
fi
