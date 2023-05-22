#!/usr/bin/env perl

# Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Reads an assembly_summary.txt file, which indicates taxids and FTP paths for
# genome/protein data.  Performs the download of the complete genomes from
# that file, decompresses, and explicitly assigns taxonomy as needed.

use strict;
use warnings;
use File::Basename;
use Getopt::Std;
use Net::FTP;
use List::Util qw/max/;

my $PROG = basename $0;
my $SERVER = "ftp.ncbi.nlm.nih.gov";
my $SERVER_PATH = "/genomes";
my $FTP_USER = "anonymous";
my $FTP_PASS = "kraken2download";

my $qm_server = quotemeta $SERVER;
my $qm_server_path = quotemeta $SERVER_PATH;

my $is_protein = $ENV{"KRAKEN2_PROTEIN_DB"};
my $use_ftp = $ENV{"KRAKEN2_USE_FTP"};

my $suffix = $is_protein ? "_protein.faa.gz" : "_genomic.fna.gz";

# Manifest hash maps filenames (keys) to taxids (values)
my %manifest;
while (<>) {
  next if /^#/;
  chomp;
  my @fields = split /\t/;
  my ($taxid, $asm_level, $ftp_path) = @fields[5, 11, 19];
  # Possible TODO - make the list here configurable by user-supplied flags
  next unless grep {$asm_level eq $_} ("Complete Genome", "Chromosome");
  next if $ftp_path eq "na";  # Skip if no provided path

  my $full_path = $ftp_path . "/" . basename($ftp_path) . $suffix;
  # strip off server/leading dir name to allow --files-from= to work w/ rsync
  # also allows filenames to just start with "all/", which is nice
  if (! ($full_path =~ s#^(?:ftp|https)://${qm_server}${qm_server_path}/##)) {
    die "$PROG: unexpected FTP path (new server?) for $ftp_path\n";
  }
  $manifest{$full_path} = $taxid;
}

open MANIFEST, ">", "manifest.txt"
  or die "$PROG: can't write manifest: $!\n";
print MANIFEST "$_\n" for keys %manifest;
close MANIFEST;

if ($is_protein && ! $use_ftp) {
  print STDERR "Step 0/2: performing rsync dry run (only protein d/l requires this)...\n";
  # Protein files aren't always present, so we have to do this two-rsync run hack
  # First, do a dry run to find non-existent files, then delete them from the
  # manifest; after this, execution can proceed as usual.
  system("rsync --dry-run --no-motd --files-from=manifest.txt rsync://${SERVER}${SERVER_PATH} . 2> rsync.err");
  open ERR_FILE, "<", "rsync.err"
    or die "$PROG: can't read rsync.err file: $!\n";
  while (<ERR_FILE>) {
    chomp;
    # I really doubt this will work across every version of rsync. :(
    if (/failed: No such file or directory/ && /^rsync: link_stat "\/([^"]+)"/) {
      delete $manifest{$1};
    }
  }
  close ERR_FILE;
  print STDERR "Rsync dry run complete, removing any non-existent files from manifest.\n";

  # Rewrite manifest
  open MANIFEST, ">", "manifest.txt"
    or die "$PROG: can't write manifest: $!\n";
  print MANIFEST "$_\n" for keys %manifest;
  close MANIFEST;
}

sub ftp_connection {
    my $ftp = Net::FTP->new($SERVER, Passive => 1)
        or die "$PROG: FTP connection error: $@\n";
    $ftp->login($FTP_USER, $FTP_PASS)
        or die "$PROG: FTP login error: " . $ftp->message() . "\n";
    $ftp->binary()
        or die "$PROG: FTP binary mode error: " . $ftp->message() . "\n";
    $ftp->cwd($SERVER_PATH)
        or die "$PROG: FTP CD error: " . $ftp->message() . "\n";
    return $ftp;
}

if ($use_ftp) {
  print STDERR "Step 1/2: Performing ftp file transfer of requested files\n";
  open MANIFEST, "<", "manifest.txt"
    or die "$PROG: can't open manifest: $!\n";
  mkdir "all" or die "$PROG: can't create 'all' directory: $!\n";
  chdir "all" or die "$PROG: can't chdir into 'all' directory: $!\n";
  while (<MANIFEST>) {
    chomp;
    my $ftp = ftp_connection();
    my $try = 0;
    my $ntries = 5;
    my $sleepsecs = 3;
    while($try < $ntries) {
        $try++;
        last if $ftp->get($_);
        warn "$PROG: unable to download $_ on try $try of $ntries: ".$ftp->message()."\n";
        last if $try == $ntries;
        sleep $sleepsecs;
        $sleepsecs *= 3;
    }
    die "$PROG: unable to download ftp://${SERVER}${SERVER_PATH}/$_\n" if $try == $ntries;
    $ftp->quit;
  }
  close MANIFEST;
  chdir ".." or die "$PROG: can't return to correct directory: $!\n";
}
else {
  print STDERR "Step 1/2: Performing rsync file transfer of requested files\n";
  system("rsync --no-motd --files-from=manifest.txt rsync://${SERVER}${SERVER_PATH}/ .") == 0
    or die "$PROG: rsync error, exiting: $?\n";
  print STDERR "Rsync file transfer complete.\n";
}
print STDERR "Step 2/2: Assigning taxonomic IDs to sequences\n";
my $output_file = $is_protein ? "library.faa" : "library.fna";
open OUT, ">", $output_file
  or die "$PROG: can't write $output_file: $!\n";
my $projects_added = 0;
my $sequences_added = 0;
my $ch_added = 0;
my $ch = $is_protein ? "aa" : "bp";
my $max_out_chars = 0;
for my $in_filename (sort keys %manifest) {
  my $taxid = $manifest{$in_filename};
  if ($use_ftp) {  # FTP downloading doesn't create full path locally
    $in_filename = "all/" . basename($in_filename);
  }
  open IN, "gunzip -c $in_filename |" or die "$PROG: can't read $in_filename: $!\n";
  while (<IN>) {
    if (/^>/) {
      s/^>/>kraken:taxid|$taxid|/;
      $sequences_added++;
    }
    else {
      $ch_added += length($_) - 1;
    }
    print OUT;
  }
  close IN;
  unlink $in_filename;
  $projects_added++;
  my $out_line = progress_line($projects_added, scalar keys %manifest, $sequences_added, $ch_added) . "...";
  $max_out_chars = max(length($out_line), $max_out_chars);
  my $space_line = " " x $max_out_chars;
  print STDERR "\r$space_line\r$out_line" if -t STDERR;
}
close OUT;
print STDERR " done.\n" if -t STDERR;

print STDERR "All files processed, cleaning up extra sequence files...";
system("rm -rf all/") == 0
  or die "$PROG: can't clean up all/ directory: $?\n";
print STDERR " done, library complete.\n";

sub progress_line {
  my ($projs, $total_projs, $seqs, $chs) = @_;
  my $line = "Processed ";
  $line .= ($projs == $total_projs) ? "$projs" : "$projs/$total_projs";
  $line .= " project" . ($total_projs > 1 ? 's' : '') . " ";
  $line .= "($seqs sequence" . ($seqs > 1 ? 's' : '') . ", ";
  my $prefix;
  my @prefixes = qw/k M G T P E/;
  while (@prefixes && $chs >= 1000) {
    $prefix = shift @prefixes;
    $chs /= 1000;
  }
  if (defined $prefix) {
    $line .= sprintf '%.2f %s%s)', $chs, $prefix, $ch;
  }
  else {
    $line .= "$chs $ch)";
  }
  return $line;
}
