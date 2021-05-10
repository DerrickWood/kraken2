#!/usr/bin/env perl

# Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Parses RDP sequence data to create Kraken taxonomy
# and sequence ID -> taxonomy ID mapping
# Input (as <>): current_{Archaea,Bacteria}_unaligned.fa

use strict;
use warnings;
use File::Basename;

my $PROG = basename $0;

my %seqid_map;
my %seen_it;
my %child_data = ("root;no rank" => {});
LINE: while (<>) {
  next unless s/^>//;
  chomp;
  my ($seq_label, $taxo_str) = split /\t/;
  my ($seqid) = split " ", $seq_label;
  $taxo_str =~ s/^Lineage=Root;rootrank;/root;no rank;/;
  $taxo_str =~ s/;$/;no rank/;  # adjust for unclassified things
  $seqid_map{$seqid} = $taxo_str;
  next if $seen_it{$taxo_str}++;
  while ($taxo_str =~ s/(;[^;]+;[^;]+)$//) {
    my $level = $1;
    my $parent = $taxo_str;
    $child_data{$parent} ||= {};
    $child_data{$parent}->{"$taxo_str$level"}++;
    next LINE if $seen_it{$taxo_str}++;
  }
}

# Assign IDs through BFS of tree, report names/nodes info in
# NCBI format
my %id_map;
my $next_node_id = 1;
open NAMES, ">", "names.dmp" or die "$PROG: can't write names.dmp: $!\n";
open NODES, ">", "nodes.dmp" or die "$PROG: can't write nodes.dmp: $!\n";
my @bfs_queue = (["root;no rank", 1]);
while (@bfs_queue) {
  my $arg_ref = shift @bfs_queue;
  my ($node, $parent_id) = @$arg_ref;
  if ($node !~ /([^;]+);([^;]+)$/) {
    die "$PROG: BFS processing encountered formatting error, \"$node\"\n";
  }
  my ($display_name, $rank) = ($1, $2);
  $rank = "superkingdom" if $rank eq "domain";  # conform w/ NCBI taxonomy
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
for my $seqid (sort keys %seqid_map) {
  my $taxid = $id_map{ $seqid_map{$seqid} };
  print SEQID_TAXID_MAP "$seqid\t$taxid\n";
}
close SEQID_TAXID_MAP;
