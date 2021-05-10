#!/bin/bash

# Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Build the standard Kraken database
# Designed to be called by kraken_build

set -u  # Protect against uninitialized vars.
set -e  # Stop on error
set -o pipefail  # Stop on failures in non-final pipeline commands

protein_flag=""
if [ -n "$KRAKEN2_PROTEIN_DB" ]; then
  protein_flag="--protein"
fi

masking_flag=""
if [ -z "$KRAKEN2_MASK_LC" ]; then
  masking_flag="--no-mask"
fi

ftp_flag=""
if [ -n "$KRAKEN2_USE_FTP" ]; then
  ftp_flag="--use-ftp"
fi

kraken2-build --db $KRAKEN2_DB_NAME --download-taxonomy $masking_flag $protein_flag $ftp_flag
kraken2-build --db $KRAKEN2_DB_NAME --download-library archaea $masking_flag $protein_flag $ftp_flag
kraken2-build --db $KRAKEN2_DB_NAME --download-library bacteria $masking_flag $protein_flag $ftp_flag
kraken2-build --db $KRAKEN2_DB_NAME --download-library viral $masking_flag $protein_flag $ftp_flag
kraken2-build --db $KRAKEN2_DB_NAME --download-library plasmid $masking_flag $protein_flag $ftp_flag
kraken2-build --db $KRAKEN2_DB_NAME --download-library human --no-mask $protein_flag $ftp_flag
if [ -z "$KRAKEN2_PROTEIN_DB" ]; then
  kraken2-build --db $KRAKEN2_DB_NAME --download-library UniVec_Core $masking_flag $ftp_flag
fi
kraken2-build --db $KRAKEN2_DB_NAME --build --threads $KRAKEN2_THREAD_CT \
              --minimizer-len $KRAKEN2_MINIMIZER_LEN \
              --kmer-len $KRAKEN2_KMER_LEN \
              --minimizer-spaces $KRAKEN2_MINIMIZER_SPACES \
              $protein_flag
