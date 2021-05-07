#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $PROG = basename $0;

my $random_seed = 42;
my $output_format = "reads#.fq";
my $num_frags = 100;
my $frag_dist_params = "300,30";
my $read_length = 100;
my $error_rate = 0.01;

GetOptions(
  "random-seed:i" => \$random_seed,
  "output-format:s" => \$output_format,
  "num-frags:i" => \$num_frags,
  "frag-dist-params:s" => \$frag_dist_params,
  "read-length:i" => \$read_length,
  "error_rate:f" => \$error_rate,
) or die "$PROG: option parsing error\n";

my ($frag_dist_mu, $frag_dist_sd) = split /,/, $frag_dist_params;
die "$PROG: illegal fragment length distribution $frag_dist_params\n" if $frag_dist_mu <= 0 || $frag_dist_sd * 6 > $frag_dist_mu;

if ($random_seed >= 0) {
  srand($random_seed);
}

my $paired = $output_format =~ /#/;
my ($fh1, $fh2);
if ($paired) {
  my $f1 = $output_format;
  $f1 =~ s/#/_1/;
  open $fh1, ">", $f1 or die "$PROG: can't write $f1: $!\n";
  my $f2 = $output_format;
  $f2 =~ s/#/_2/;
  open $fh2, ">", $f2 or die "$PROG: can't write $f2: $!\n";
}
else {
  open $fh1, ">", $output_format
    or die "$PROG: can't write $output_format: $!\n";
}

my %references;
my ($seqid, $seq);
while (<>) {
  if (/^>(\S+)/) {
    $references{$seqid} = $seq if defined $seqid;
    $seqid = $1;
    $seq = "";
  }
  else {
    chomp;
    $seq .= $_;
  }
}
$references{$seqid} = $seq if defined $seqid;
my $total_reference_length = 0;
$total_reference_length += length($_) for values %references;

for my $fragment_id (1 .. $num_frags) {
  my $ref_id = random_reference();
  my $ref = $references{$ref_id};
  my $fraglen;
  do { $fraglen = int(rand_normal() * $frag_dist_sd + $frag_dist_mu) } until $fraglen <= length($ref);
  my $startpos = int rand(length($ref) - $fraglen + 1);
  my $fragment = substr($ref, $startpos, $fraglen);
  my $seq1 = substr($fragment, 0, $read_length);
  my $seq2 = revcom(substr($fragment, -$read_length, $read_length));
  my ($read1, $qual1) = errorize($seq1);
  my ($read2, $qual2) = errorize($seq2);
  print $fh1 "\@FRAG.${fragment_id}/1 $ref_id:@{[$startpos + 1]}\n$read1\n+\n$qual1\n";
  print $fh2 "\@FRAG.${fragment_id}/2 $ref_id:@{[$startpos + 1]}\n$read2\n+\n$qual2\n" if $paired;
}

sub errorize {
  my $seq = shift;
  my $quals = "J" x length($seq);
  for my $i (0 .. length($seq) - 1) {
    if (rand() < $error_rate) {
      my @ch = grep { substr($seq, $i, 1) ne $_ } qw/A C G T/;
      my $edit = $ch[rand @ch];
      substr($seq, $i, 1) = $edit;
      substr($quals, $i, 1) = '#';
    }
  }
  while (length($seq) < $read_length) {
    $seq .= "N";
    $quals .= "#";
  }
  return ($seq, $quals);
}

sub random_reference {
  # Try to do random selection weighted by sequence length
  my $r = rand();
  my $p = 0;
  while (my ($id, $seq) = each %references) {
    $p += length($seq) / $total_reference_length;
    if ($r < $p) {
      return $id;
    }
  }
  # Just return random ID as fallback
  my @k = keys %references;
  return $k[rand @k];
}

{ # Return random numbers with standard Gaussian distribution using Box-Muller
  my ($r, $th);
  sub rand_normal {
    my $z;
    if (defined $r) {
      $z = $r * sin($th);
      undef($r);
    }
    else {
      $r = sqrt(-2 * log(rand()));
      $th = 2 * 3.1415926 * rand();
      $z = $r * cos($th);
    }
    return $z;
  }
}

sub revcom {
  local $_ = shift;
  tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
  return scalar reverse;
}
