#!/usr/bin/env perl

# Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Create a file in a specified directory, then copy an
# existing file's contents into the new file.  Write name of
# new file to standard output.
#
# This exists because the mktemp program doesn't act consistently across
# operating systems/distros/versions.

use strict;
use warnings;
use File::Basename;
use File::Temp 'tempfile';
use Getopt::Std;

my $PROG = basename $0;
getopts('d:t:s:', \my %opts) or usage();
$opts{$_} or usage() for qw/d t s/;  # all switches mandatory
my ($directory, $template, $suffix) = @opts{qw/d t s/};
die "$PROG: '$directory' not a directory!\n" unless -d $directory;
die "$PROG: must specify a single filename\n" unless @ARGV == 1;

$suffix =~ s/^\.//;
my $old_filename = shift @ARGV;
open FILE, "<", $old_filename
  or die "$PROG: can't read $old_filename: $!\n";

my ($fh, $new_filename) = tempfile($template, DIR => $directory,
                                   UNLINK => 0, SUFFIX => ".$suffix");
# copy loop
while (<FILE>) {
  print {$fh} $_;
}
close FILE;
close $fh;

print "$new_filename\n";

sub usage {
  die "$PROG: <-d directory> <-t template> <-s suffix> <filename>\n";
}
