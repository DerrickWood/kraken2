#!/bin/bash
# Regenerates the frozen read fixtures in tests/data from data/*.fa.
#
# The reads are committed rather than simulated at test time because
# data/simulator.pl walks a Perl hash to pick references, and Perl hash
# order varies between runs and Perl versions.  Pinning the hash seed makes
# this script reproducible on one Perl build, which is all that is needed to
# regenerate the fixtures deliberately.  Changing the reads invalidates every
# golden file, so run tests/update_golden.sh afterwards.

set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
cd "$here"

PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0 perl "$repo/data/simulator.pl" \
  --random-seed 42 --num-frags 400 --output-format sim#.fq "$repo"/data/*.fa

# Append chimeric pairs whose mates come from two different taxa.  Each mate
# is an exact 100 bp fragment whose k-mers all map to its own leaf taxon, so
# the two root-to-leaf paths score the same and the call must be the LCA of
# the two taxa.  These exercise the tie-breaking branch of ResolveTree, which
# ordinary simulated reads almost never reach.  The candidates come from a
# seeded draw and KEEP lists the ones verified to be exact ties.
perl -e '
  use strict;
  my %keep = map { $_ => 1 } qw(2 3 4 10 15 18 19 22 29 33 49 51 53 56 62 65 80 82 88 100);
  my $repo = shift;
  my (%parent, %seq, %tax, $id);
  open my $n, "<", "$repo/data/nodes.dmp" or die;
  while (<$n>) { my @f = split /\t\|\t/; $parent{$f[0]} = $f[1] }
  for my $f (sort glob "$repo/data/*.fa") {
    open my $fh, "<", $f or die;
    while (<$fh>) {
      chomp;
      if (/^>(\S+)/) { $id = $1; ($tax{$id}) = $id =~ /taxid\|(\d+)\|/ }
      else { $seq{$id} .= $_ }
    }
  }
  sub lineage { my @l = (shift); push @l, $parent{$l[-1]} while $l[-1] != 1; @l }
  sub lca { my %a = map { $_ => 1 } lineage(shift); for (lineage(shift)) { return $_ if $a{$_} } }
  my @ids = sort keys %seq;
  srand(11);
  open my $f1, ">>", "sim_1.fq" or die;
  open my $f2, ">>", "sim_2.fq" or die;
  for my $i (1 .. 100) {
    my ($a, $b) = ($ids[int rand @ids], $ids[int rand @ids]);
    next if $tax{$a} == $tax{$b};
    my $pa = int rand(length($seq{$a}) - 100);
    my $pb = int rand(length($seq{$b}) - 100);
    next unless $keep{$i};
    my $sa = substr($seq{$a}, $pa, 100);
    my $sb = reverse substr($seq{$b}, $pb, 100);
    $sb =~ tr/ACGT/TGCA/;
    my $t = lca($tax{$a}, $tax{$b});
    my $q = "J" x 100;
    print $f1 "\@CHIMERA.$i/1 kraken:taxid|$t|chimera $a:$pa\n$sa\n+\n$q\n";
    print $f2 "\@CHIMERA.$i/2 kraken:taxid|$t|chimera $b:$pb\n$sb\n+\n$q\n";
  }' "$repo"

# Append pairs of random sequence that should stay unclassified.  Their
# names start with RANDOM so the accuracy check treats them as taxid 0.
perl -e '
  srand(7);
  open my $f1, ">>", "sim_1.fq" or die;
  open my $f2, ">>", "sim_2.fq" or die;
  for my $i (1 .. 20) {
    for my $fh ($f1, $f2) {
      my $s = join "", map { (qw/A C G T/)[int rand 4] } 1 .. 100;
      my $m = $fh == $f1 ? 1 : 2;
      print $fh "\@RANDOM.$i/$m random:0\n$s\n+\n", "J" x 100, "\n";
    }
  }'
