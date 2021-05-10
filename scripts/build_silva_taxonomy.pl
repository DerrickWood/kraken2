#!/usr/bin/env perl

# Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Parses Silva taxonomy file to create Kraken taxonomy
# Input (as <>): tax_slv_ssu_nr_119.txt

use strict;
use warnings;
use File::Basename;

my $PROG = basename $0;

my %id_map = ("root" => 1);
open NAMES, ">", "names.dmp" or die "$PROG: can't write names.dmp: $!\n";
open NODES, ">", "nodes.dmp" or die "$PROG: can't write nodes.dmp: $!\n";
print NAMES "1\t|\troot\t|\t-\t|\tscientific name\t|\n";
print NODES "1\t|\t1\t|\tno rank\t|\t-\t|\n";
while (<>) {
  chomp;
  my ($taxo_str, $node_id, $rank) = split /\t/;
  $id_map{$taxo_str} = $node_id;
  if ($taxo_str =~ /^(.+;|)([^;]+);$/) {
    my $parent_name = $1;
    my $display_name = $2;
    if ($parent_name eq "") {
      $parent_name = "root";
    }
    my $parent_id = $id_map{$parent_name};
    if (! defined $parent_id) {
      die "$PROG: orphan error, line $.\n";
    }
    $rank = "superkingdom" if $rank eq "domain";
    print NAMES "$node_id\t|\t$display_name\t|\t-\t|\tscientific name\t|\n";
    print NODES "$node_id\t|\t$parent_id\t|\t$rank\t|\t-\t|\n";
  }
  else {
    die "$PROG: strange input, line $.\n";
  }
}
close NAMES;
close NODES;
