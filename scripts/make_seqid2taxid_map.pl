#!/usr/bin/env perl

# Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Reads multi-FASTA input and examines each sequence header.  In quiet mode
# (-q, used as an initial check on file validity), headers are OK if a
# taxonomy ID is found (as either the entire sequence ID or as part of a
# "kraken:taxid" token), or if something looking like a GI or accession
# number is found.  In normal mode, the taxonomy ID will be looked up (if
# not explicitly specified in the sequence ID) and reported if it can be
# found.  Output is tab-delimited lines, with sequence IDs in first
# column and taxonomy IDs in second.

# Sequence IDs with a kraken:taxid token will use that to assign taxonomy
# ID, e.g.:
# >gi|32499|ref|NC_021949.2|kraken:taxid|562|
#
# Sequence IDs that are completely numeric are assumed to be the taxonomy
# ID for that sequence.
#
# Otherwise, an accession number is searched for; if not found, a GI
# number is searched for.  Failure to find any of the above is a fatal error.
# Without -q, a comma-separated file list specified by -A (for both accession
# numbers and GI numbers) is examined; failure to find a
# taxonomy ID that maps to a provided accession/GI number is non-fatal and
# will emit a warning.
#
# With -q, does not print any output, and will die w/ nonzero exit instead
# of warning when unable to find a taxid, accession #, or GI #.

use strict;
use warnings;
use File::Basename;
use Getopt::Std;

my $PROG = basename $0;
getopts('qA:L:', \my %opts);

if ($opts{'q'} && (defined $opts{"A"} || defined $opts{"L"})) {
  die "$PROG: -q doesn't allow for -A/-L\n";
}

my %target_lists;

while (<>) {
  next unless /^>(\S+)/;
  my $seq_id = $1;
  my $output;
  # Note all regexes here use ?: to avoid capturing the ^ or | character
  if ($seq_id =~ /(?:^|\|)kraken:taxid\|(\d+)/) {
    $output = "$seq_id\t$1\n";
  }
  elsif ($seq_id =~ /^\d+$/) {
    $output = "$seq_id\t$seq_id\n";
  }
  # Accession regex is first, only puts number (not version) into $1
  elsif ($seq_id =~ /(?:^|\|)([A-Z]+_?[A-Z0-9]+)(?:\||\b|\.)/ ||
         $seq_id =~ /(?:^|\|)gi\|(\d+)/)
  {
    if (! $opts{"q"}) {
      $target_lists{$1} ||= [];
      push @{ $target_lists{$1} }, $seq_id;
    }
  }
  else {
    die "$PROG: unable to determine taxonomy ID for sequence $seq_id\n";
  }

  if (defined $output && ! $opts{"q"}) {
    print $output;
  }
}

if ($opts{"q"}) {
  if (keys %target_lists) {
    print "$PROG: requires external map\n";
  }
  exit 0;
}
exit 0 if ! keys %target_lists;
if (! (defined $opts{"A"} || defined $opts{"L"})) {
  die "$PROG: found sequence ID without explicit taxonomy ID, but no map used\n";
}

# Remove targets where we've already handled the mapping
if (defined $opts{"L"}) {
  my $library_map_file = $opts{"L"};
  open FILE, "<", $library_map_file
    or die "$PROG: can't open $library_map_file: $!\n";
  while (<FILE>) {
    chomp;
    my ($seqid, $taxid) = split /\t/;
    if (exists $target_lists{$seqid}) {
      print "$seqid\t$taxid\n";
      delete $target_lists{$seqid};
    }
  }
  close FILE;
}

exit 0 if ! keys %target_lists;

my @accession_map_files = split /,/, $opts{"A"};
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
    }
    if ($gi ne "na" && $target_lists{$gi}) {
      my @list = @{ $target_lists{$gi} };
      delete $target_lists{$gi};
      for my $seqid (@list) {
        print "$seqid\t$taxid\n";
      }
    }
  }
  close FILE;
}
