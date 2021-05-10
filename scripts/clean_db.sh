#!/bin/bash

# Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Removes intermediate files from a database directory,
# such as reference library FASTA files and taxonomy data from NCBI.

set -u  # Protect against uninitialized vars.
set -e  # Stop on error

cd $KRAKEN2_DB_NAME
previous_usage=$(du -sh | cut -f1)
1>&2 echo "Database disk usage: $previous_usage"
rm -rf library/ taxonomy/ seqid2taxid.map
current_usage=$(du -sh | cut -f1)
1>&2 echo "After cleaning, database uses $current_usage"
