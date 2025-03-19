#!/bin/sh

# Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

set -e

SCRIPT="$(realpath "$0")"
ROOT=$(dirname "$SCRIPT")
VERSION="$(cat "$ROOT/VERSION")"

cd "$ROOT"

if [ -z "$1" ] || [ -n "$2" ]
then
  echo "Usage: $(basename "$SCRIPT") KRAKEN2_DIR"
  exit 64
fi

if [ "$1" = "KRAKEN2_DIR" ]
then
  echo "Please replace \"KRAKEN2_DIR\" with the name of the directory"
  echo "that you want to install Kraken 2 in."
  exit 1
fi

# Perl cmd used to canonicalize dirname - "readlink -f" doesn't work
# on OS X.
export KRAKEN2_DIR
KRAKEN2_DIR=$(perl -MCwd=abs_path -le 'print abs_path(shift)' "$1")

mkdir -p "$KRAKEN2_DIR"
make -C src install
for file in scripts/*
do
  destination_file="$KRAKEN2_DIR/$(basename "$file")"
  perl -pl -e 'BEGIN { while (@ARGV) { $_ = shift; ($k,$v) = split /=/, $_, 2; $H{$k} = $v } }'\
           -e 's/#####=(\w+)=#####/$H{$1}/g' \
           "KRAKEN2_DIR=$KRAKEN2_DIR" "VERSION=$VERSION" \
           < "$file" > "$destination_file"
  if [ -x "$file" ]
  then
    chmod +x "$destination_file"
  fi
done

echo
echo "Kraken 2 installation complete."
echo
echo "To make things easier for you, you may want to copy/symlink the following"
echo "files into a directory in your PATH:"
for file in "$KRAKEN2_DIR"/kraken2*
do
  if [ -x "$file" ]
  then
    echo "  $file"
  fi
done
