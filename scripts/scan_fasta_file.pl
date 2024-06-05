#!/usr/bin/env perl

# Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Reads multi-FASTA input and examines each sequence header.  Headers are
# OK if a taxonomy ID is found (as either the entire sequence ID or as part
# of a # "kraken:taxid" token), or if something looking like an accession
# number is found.  Not "OK" headers will are fatal errors unless "--lenient"
# is used.
#
# Each sequence header results in a line with three tab-separated values;
# the first indicating whether third column is the taxonomy ID ("TAXID") or
# an accession number ("ACCNUM") for the sequence ID listed in the second
# column.

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
require "$ENV{KRAKEN2_DIR}/kraken2lib.pm";

my $PROG = basename $0;

my $lenient = 0;
GetOptions("lenient" => \$lenient)
  or die "Usage: $PROG [--lenient] <fasta filename(s)>\n";

while (<>) {
  next unless /^>/;
  # while (/.../g) needed because non-redundant DBs sometimes have multiple
  #   sequence IDs in the header; extra sequence IDs are prefixed by
  #   '\x01' characters (if downloaded in FASTA format from NCBI FTP directly).
  while (/(?:^>|\x01)(\S+)/g) {
    my $seqid = $1;
    my $taxid = kraken2lib::check_seqid($seqid);
    if (! defined $taxid) {
      next if $lenient;
      die "$PROG: unable to determine taxonomy ID for sequence $seqid\n";
    }
    # "$taxid" may actually be an accession number
    if ($taxid =~ /^\d+$/) {
      print "TAXID\t$seqid\t$taxid\n";
    }
    else {
      print "ACCNUM\t$seqid\t$taxid\n";
    }
  }
}
