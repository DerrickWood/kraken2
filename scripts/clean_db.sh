#!/bin/bash

# Copyright 2013-2018, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Removes intermediate files from a database directory,
# such as reference library FASTA files and taxonomy data from NCBI.

set -u  # Protect against uninitialized vars.
set -e  # Stop on error

cd $KRAKEN2_DB_NAME
rm -rf library/ taxonomy/ seqid2taxid.map
