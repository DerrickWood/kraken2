#!/usr/bin/env perl

# Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Looks up accession numbers and reports associated taxonomy IDs
#
# Input is (a) 1 2-column TSV file w/ sequence IDs and accession numbers,
# and (b) a list of accession2taxid files from NCBI.
# Output is tab-delimited lines, with sequence IDs in first
# column and taxonomy IDs in second.

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

my $PROG = basename $0;

my $lookup_list_file = shift @ARGV;

my %target_lists;
open FILE, "<", $lookup_list_file
  or die "$PROG: can't open $lookup_list_file: $!\n";
while (<FILE>) {
  chomp;
  my ($seqid, $accnum) = split /\t/;
  $target_lists{$accnum} ||= [];
  push @{ $target_lists{$accnum} }, $seqid;
}
close FILE;

my $initial_target_count = scalar keys %target_lists;

my @accession_map_files = @ARGV;
for my $file (@accession_map_files) {
  open FILE, "<", $file
    or die "$PROG: can't open $file: $!\n";
  scalar(<FILE>);  # discard header line
  while (<FILE>) {
    chomp;
    my ($accession, $with_version, $taxid, $gi) = split /\t/;
    if ($target_lists{$accession}) {
      my @list = @{ $target_lists{$accession} };
      delete $target_lists{$accession};
      for my $seqid (@list) {
        print "$seqid\t$taxid\n";
      }
      last if ! %target_lists;
    }
  }
  close FILE;
  last if ! %target_lists;
}

if (%target_lists) {
  warn "$PROG: @{[scalar keys %target_lists]}/$initial_target_count accession numbers remain unmapped, see unmapped.txt in DB directory\n";
  open LIST, ">", "unmapped.txt"
    or die "$PROG: can't write unmapped.txt: $!\n";
  for my $target (keys %target_lists) {
    print LIST "$target\n";
  }
  close LIST;
}
