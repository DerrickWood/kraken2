package kraken2lib;

# Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Common subroutines for other Kraken scripts

use strict;
use warnings;

# Input: the argument for a --db option (possibly undefined)
# Returns: the DB to use, taking KRAKEN2_DEFAULT_DB and KRAKEN2_DB_PATH
#   into account.
sub find_db {
  my $supplied_db_prefix = shift;
  my $db_prefix;
  if (! defined $supplied_db_prefix) {
    if (! exists $ENV{"KRAKEN2_DEFAULT_DB"}) {
      die "Must specify DB with either --db or \$KRAKEN2_DEFAULT_DB\n";
    }
    $supplied_db_prefix = $ENV{"KRAKEN2_DEFAULT_DB"};
  }
  my @db_path = (".");
  if (exists $ENV{"KRAKEN2_DB_PATH"}) {
    my $path_str = $ENV{"KRAKEN2_DB_PATH"};
    # Allow zero-length path to be current dir
    $path_str =~ s/^:/.:/;
    $path_str =~ s/:$/:./;
    $path_str =~ s/::/:.:/;

    @db_path = split /:/, $path_str;
  }

  # Use supplied DB if abs. or rel. path is given
  if ($supplied_db_prefix =~ m|/|) {
    $db_prefix = $supplied_db_prefix;
  }
  else {
    # Check all dirs in KRAKEN2_DB_PATH
    for my $dir (@db_path) {
      my $checked_db = "$dir/$supplied_db_prefix";
      if (-e $checked_db && -d _) {
        $db_prefix = $checked_db;
        last;
      }
    }
    if (! defined $db_prefix) {
      my $printed_path = exists $ENV{"KRAKEN2_DB_PATH"} ? qq|"$ENV{'KRAKEN2_DB_PATH'}"| : "undefined";
      die "unable to find $supplied_db_prefix in \$KRAKEN2_DB_PATH ($printed_path)\n";
    }
  }

  for my $file (qw/taxo.k2d hash.k2d opts.k2d/) {
    if (! -e "$db_prefix/$file") {
      die "database (\"$db_prefix\") does not contain necessary file $file\n";
    }
  }

  return $db_prefix;
}

# Input: a FASTA sequence ID
# Output: either (a) a taxonomy ID number found in the sequence ID,
#   (b) an NCBI accession number found in the sequence ID, or undef
sub check_seqid {
  my $seqid = shift;
  my $taxid = undef;
  # Note all regexes here use ?: to avoid capturing the ^ or | character
  if ($seqid =~ /(?:^|\|)kraken:taxid\|(\d+)/) {
    $taxid = $1;  # OK, has explicit taxid
  }
  elsif ($seqid =~ /^(\d+)$/) {
    $taxid = $1;  # OK, has explicit taxid (w/o token)
  }
  # Accession number check
  elsif ($seqid =~ /(?:^|\|)         # Begins seqid or immediately follows pipe
                     ([A-Z]+         # Starts with one or more UC alphas
                        _?           # Might have an underscore next
                        [A-Z0-9]+)   # Ends with UC alphas or digits
                     (?:\||\b|\.)/x  # Followed by pipe, word boundary, or period
        )
  {
    $taxid = $1;  # A bit misleading - failure to pass /^\d+$/ means this is
                  # OK, but requires accession -> taxid mapping
  }
  return $taxid;
}

1;
