#!/usr/bin/env perl

# Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Parses Greengenes taxonomy file to create Kraken taxonomy
# and sequence ID -> taxonomy ID mapping
# Input (as <>): gg_13_5_taxonomy.txt

use strict;
use warnings;
use File::Basename;

my $PROG = basename $0;

my %RANK_CODES = (
  k => "superkingdom",
  p => "phylum",
  c => "class",
  o => "order",
  f => "family",
  g => "genus",
  s => "species"
);

my %seqid_map;
my %seen_it;
my %child_data = ("root" => {});
LINE: while (<>) {
  chomp;
  my ($seqid, $taxo_str) = split /\t/;
  $taxo_str =~ s/(; [a-z]__)+$//;  # Remove empty data
  $seqid_map{$seqid} = $taxo_str;
  next if $seen_it{$taxo_str}++;
  while ($taxo_str =~ s/(; [a-z]__[^;]+$)//) {
    my $level = $1;
    my $parent = $taxo_str;
    $child_data{$parent} ||= {};
    $child_data{$parent}->{"$taxo_str$level"}++;
    next LINE if $seen_it{$taxo_str}++;
  }
  $child_data{"root"}->{$taxo_str}++;
}

# Assign IDs through BFS of tree, report names/nodes info in
# NCBI format
my %id_map;
my $next_node_id = 1;
open NAMES, ">", "names.dmp" or die "$PROG: can't write names.dmp: $!\n";
open NODES, ">", "nodes.dmp" or die "$PROG: can't write nodes.dmp: $!\n";
my @bfs_queue = (["root", 1]);
while (@bfs_queue) {
  my $arg_ref = shift @bfs_queue;
  my ($node, $parent_id) = @$arg_ref;
  my $display_name = $node;
  my $rank;
  # special handling for species
  if ($node =~ /g__([^;]+); s__([^;]+)$/) {
    my ($genus, $species) = ($1, $2);
    $rank = "species";
    if ($species =~ / endosymbiont /) {
      $display_name = $species;
    }
    else {
      $display_name = "$genus $species";
    }
  }
  elsif ($node =~ /([a-z])__([^;]+)$/) {
    $rank = $RANK_CODES{$1};
    $display_name = $2;
  }
  $rank ||= "no rank";
  my $node_id = $next_node_id++;
  $id_map{$node} = $node_id;
  print NAMES "$node_id\t|\t$display_name\t|\t-\t|\tscientific name\t|\n";
  print NODES "$node_id\t|\t$parent_id\t|\t$rank\t|\t-\t|\n";
  my @children = sort keys %{ $child_data{$node} };
  push @bfs_queue, [$_, $node_id] for @children;
}
close NAMES;
close NODES;

open SEQID_TAXID_MAP, ">", "seqid2taxid.map" or die "$PROG: can't write seqid2taxid.map: $!\n";
for my $seqid (sort { $a <=> $b } keys %seqid_map) {
  my $taxid = $id_map{ $seqid_map{$seqid} };
  print SEQID_TAXID_MAP "$seqid\t$taxid\n";
}
close SEQID_TAXID_MAP;
