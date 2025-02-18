#!/usr/bin/env python3

import argparse
import bz2
import concurrent.futures
import collections
import ftplib
import glob
import gzip
import hashlib
import http.client
import inspect
import io
import itertools
import json
import logging
import lzma
import math
import os
import pty
import random
import re
import shutil
import signal
import subprocess
import sys
import tarfile
import tempfile
import threading
import urllib
import urllib.error
import urllib.parse
import urllib.request
import zipfile
import zlib

LOG = None
SCRIPT_PATHNAME = None

NCBI_REST_API = "api.ncbi.nlm.nih.gov"
NCBI_SERVER = "ftp.ncbi.nlm.nih.gov"
GREENGENES_SERVER = "greengenes.microbio.me"
SILVA_SERVER = "ftp.arb-silva.de"
GTDB_SERVER = "data.gtdb.ecogenomic.org"

WRAPPER_ARGS_TO_BIN_ARGS = {
    "block_size": "-B",
    "classified_out": "-C",
    "confidence": "-T",
    "fast_build": "-F",
    "interleaved": "-S",
    "kmer_len": "-k",
    "max_db_size": "-M",
    "memory_mapping": "-M",
    "minimizer_len": "-l",
    "minimum_bits_for_taxid": "-r",
    "minimum_base_quality": "-Q",
    "minimum_hit_groups": "-g",
    "output": "-O",
    "paired": "-P",
    "protein": "-X",
    "quick": "-q",
    "report": "-R",
    "report_minimizer_data": "-K",
    "report_zero_counts": "-z",
    "skip_counts": "-s",
    "sub_block_size": "-b",
    "threads": "-p",
    "unclassified_out": "-U",
    "use_mpa_style": "-m",
    "use_names": "-n",
}


class FTP:
    def __init__(self, server):
        self.ftp = ftplib.FTP(server, timeout=600)
        self.ftp.login()
        self.ftp.sendcmd("TYPE I")
        self.pwd = "/"
        self.server = server

    def _progress_bar(self, f, remote_size):
        pb = ProgressBar(remote_size, f.tell())

        def inner(block):
            nonlocal f, remote_size, pb
            written = 0
            while written < len(block):
                written += f.write(block[written:])
            size_on_disk = f.tell()
            pb.progress(size_on_disk)
            LOG.debug(
                "{:s} {: >10s}\r".format(
                    pb.get_bar(), format_bytes(size_on_disk)
                )
            )

        return inner

    def download(self, remote_dir, filepaths):
        if isinstance(filepaths, str):
            filepaths = [filepaths]
        number_of_files = len(filepaths)
        self.cwd(remote_dir)
        for index, filepath in enumerate(filepaths):
            mode = "ab"
            local_size = 0
            remote_size = self.size(filepath)
            if os.path.exists(filepath):
                local_size = os.stat(filepath).st_size
            else:
                if os.path.basename(filepath) != filepath:
                    os.makedirs(os.path.dirname(filepath), exist_ok=True)
            if local_size == remote_size:
                LOG.info(
                    "Already downloaded {:s}\n".format(get_abs_path(filepath))
                )
                continue
            if local_size > remote_size:
                mode = "wb"
            url_components = urllib.parse.SplitResult(
                "ftp", self.server, os.path.join(remote_dir, filepath), "", ""
            )
            url = urllib.parse.urlunsplit(url_components)
            if number_of_files == 1:
                LOG.info("Downloading {:s}\n".format(url))
            else:
                LOG.info(
                    "[{:d}/{:d}] Downloading {:s}\n".format(
                        index + 1, number_of_files, url
                    )
                )
            with open(filepath, mode) as f:
                while True:
                    try:
                        cb = self._progress_bar(f, remote_size)
                        self.ftp.retrbinary(
                            "RETR " + filepath, cb, rest=f.tell()
                        )
                        break
                    except KeyboardInterrupt:
                        f.flush()
                        self.close()
                        sys.exit(1)
                    except ftplib.all_errors:
                        f.flush()
                        self.reconnect()
                        self.cwd(remote_dir)
                        continue
            absolute_path = get_abs_path(filepath)
            local_filename, local_dirname = os.path.basename(
                absolute_path
            ), os.path.dirname(absolute_path)
            clear_console_line()
            LOG.info(
                "Saved {:s} to {:s}\n".format(local_filename, local_dirname)
            )

    def cwd(self, remote_pathname):
        self.ftp.cwd(remote_pathname)
        self.pwd = remote_pathname

    def size(self, filepath):
        size = 0
        while True:
            try:
                size = self.ftp.size(filepath)
                break
            except ftplib.error_temp:
                self.reconnect()
                continue
        return size

    def exists(self, filepath):
        while True:
            try:
                self.size(filepath)
                break
            except ftplib.error_perm as e:
                if e.args[0].find("No such file or directory"):
                    return False
                raise
        return True

    def connect(self, server):
        self.ftp = ftplib.FTP(server)
        self.ftp.login()
        self.ftp.sendcmd("TYPE I")

    def reconnect(self):
        host = self.ftp.host
        self.ftp.close()
        self.connect(host)
        self.ftp.cwd(self.pwd)

    def host(self):
        return self.ftp.host

    def close(self):
        self.ftp.quit()


class ProgressBar:
    def __init__(self, stop, current=0, width=30):
        self.stop = stop
        self.width = width
        self.current = current
        self.bar = list("-" * self.width)
        self.step = stop / self.width
        self.last_index = self._calculate_index()
        if self.current > 0:
            self.progress()

    def progress(self, amount=0, relative=False):
        if relative:
            self.current += amount
        else:
            self.current = amount
        if self.current > self.stop:
            self.current = self.stop
        index = self._calculate_index()
        for i in range(self.last_index, index):
            if i == 0:
                self.bar[i] = ">"
            else:
                self.bar[i - 1], self.bar[i] = "=", ">"
        self.last_index = index

    def get_bar(self):
        percentage = int(self.current / self.stop * 100)
        return "{:3d}% {:s}".format(percentage, "[" + "".join(self.bar) + "]")

    def _calculate_index(self):
        return math.floor(self.current / self.step)


class NCBI_URI_Builder:
    def __init__(self, *path_components):
        path_components = list(path_components)
        for i, component in enumerate(path_components):
            if isinstance(component, list):
                component = ",".join(component)
            path_components[i] = urllib.parse.quote(component)
        self.filters = {}
        self.path = "/datasets/v2alpha/genome/{}".format("/".join(path_components))

    def assembly_source(self, source=None):
        if source:
            self.filters["filters.assembly_source"] = urllib.parse.quote(source)
        return self

    def assembly_levels(self, levels):
        if levels:
            self.filters["filters.assembly_level"] = levels
        return self

    def assembly_version(self, version=None):
        if version:
            self.filters["filters.assembly_version"] = version
        return self

    def exclude_paired_reports(self, exclude_pairs=False):
        if exclude_pairs:
            self.filters["filters.exclude_paired_reports"] = "true"
        return self

    def has_annotation(self, annotated=False):
        if annotated:
            self.filters["filters.has_annotation"] = "true"
        return self

    def search_text(self, text=None):
        if text:
            self.filters["filters.search_text"] = urllib.parse.quote(text)
        return self

    def reference_only(self, reference_only=False):
        if reference_only:
            self.filters["filters.reference_only"] = reference_only
        return self

    def page_size(self, size=None):
        if size:
            self.filters["page_size"] = size
        return self

    def page_token(self, token=None):
        if token:
            self.filters["page_token"] = token
        return self

    def include_annotation_type(self, annotation_type=None):
        if annotation_type:
            self.filters["include_annotation_type"] = annotation_type
        return self

    def set_filters_from_args(self, args):
        for k, v in vars(args).items():
            if hasattr(self, k):
                self = getattr(self, k)(v)

    def build(self):
        filters = []
        for k, v in self.filters.items():
            if isinstance(v, list):
                for value in v:
                    filters.append("{}={}".format(k, value))
            else:
                filters.append("{}={}".format(k, v))
        query = "&".join(filters)
        split = urllib.parse.SplitResult(
            scheme="",
            netloc="",
            path=self.path,
            query=query,
            fragment=""
        )

        return urllib.parse.urlunsplit(split)

    def reset(self):
        self.filters.clear()


def clear_console_line():
    LOG.debug("\33[2K\r")


def count_lines(*filenames):
    lines = 0
    for fname in filenames:
        with open(fname, "r") as f:
            for line in f:
                lines += 1
    return lines


def dwk2():
    estimate_capacity = find_kraken2_binary("estimate_capacity")
    output = subprocess.check_output(
        [estimate_capacity, "-h"], stderr=subprocess.STDOUT
    )
    for line in output.split(b"\n"):
        if line.startswith(b"Usage:"):
            return True if line.strip().endswith(b"<options>") else False
    return False


def get_binary_options(binary_pathname):
    options = []
    proc = subprocess.Popen(binary_pathname, stderr=subprocess.PIPE)
    lines = proc.stderr.readlines()
    for line in lines:
        match = re.search(rb"\s(-\w)\s", line)
        if not match:
            continue
        options.append(match.group(1).decode())
    return options


def construct_seed_template(args):
    if int(args.minimizer_len / 4) < args.minimizer_spaces:
        LOG.error(
            "Number of minimizer spaces, {}, exceeds max for"
            "minimizer length, {}; max: {}\n".format(
                args.minimizer_spaces,
                args.minimizer_len,
                int(args.minimizer_len / 4),
            )
        )
        sys.exit(1)
    return (
        "1" * (args.minimizer_len - 2 * args.minimizer_spaces)
        + "01" * args.minimizer_spaces
    )


def future_raised_exception(future):
    return future.done() and future.result() is None


def url_join(netloc, scheme="https", path="", query="", fragment=""):
    split_result = urllib.parse.SplitResult(scheme, netloc, path, query, fragment)
    return urllib.parse.urlunsplit(split_result)


def wrapper_args_to_binary_args(opts, argv, binary_args):
    for k, v in vars(opts).items():
        if k not in WRAPPER_ARGS_TO_BIN_ARGS:
            continue
        if WRAPPER_ARGS_TO_BIN_ARGS[k] not in binary_args:
            continue
        if v is False:
            continue
        if v is None:
            continue
        if v is True:
            argv.append(WRAPPER_ARGS_TO_BIN_ARGS[k])
        else:
            argv.extend([WRAPPER_ARGS_TO_BIN_ARGS[k], str(v)])


def find_kraken2_binary(name):
    # search the OS PATH
    if "PATH" in os.environ:
        for dir in os.environ["PATH"].split(":"):
            if os.path.exists(os.path.join(dir, name)):
                return os.path.join(dir, name)
    # search for binary in the same directory as wrapper
    script_parent_directory = get_parent_directory(SCRIPT_PATHNAME)
    if os.path.exists(os.path.join(script_parent_directory, name)):
        return os.path.join(script_parent_directory, name)
    # if called from within kraken2 project root, search the src dir
    project_root = get_parent_directory(script_parent_directory)
    if "src" in os.listdir(project_root) and name in os.listdir(
        os.path.join(project_root, "src")
    ):
        return os.path.join(project_root, os.path.join("src", name))
    # not found in these likely places, exit
    LOG.error("Unable to find {:s}, exiting\n".format(name))
    sys.exit(1)


def get_parent_directory(pathname):
    if len(pathname) == 0:
        return None
    pathname = os.path.abspath(pathname)
    if len(pathname) > 1 and pathname[-1] == os.path.sep:
        return os.path.dirname(pathname[:-1])
    return os.path.dirname(pathname)


def find_database(database_name):
    database_path = None
    if not os.path.isdir(database_name):
        if "KRAKEN2_DB_PATH" in os.environ:
            for directory in os.environ["KRAKEN2_DB_PATH"].split(":"):
                if os.path.exists(os.path.join(directory, database_name)):
                    database_path = os.path.join(directory, database_name)
                    break
        else:
            if database_name in os.listdir(os.getcwd()):
                database_path = database_name
    else:
        database_path = os.path.abspath(database_name)
    if database_path:
        for db_file in ["taxo.k2d", "hash.k2d", "opts.k2d"]:
            if not os.path.exists(os.path.join(database_path, db_file)):
                return None
    return database_path


def remove_files(filepaths, forked=False):
    total_size = 0

    for fname in filepaths:
        if not os.path.exists(fname):
            continue
        elif os.path.isdir(fname):
            with os.scandir(fname) as iter:
                directories = []
                for entry in iter:
                    if entry.is_dir():
                        directories.append(entry.path)
                    else:
                        filepaths.append(entry.path)
                if not forked and len(directories) >= 4:
                    total_size += remove_files_parallel(directories)
                else:
                    filepaths.extend(directories)
        else:
            LOG.info("Removing {}\n".format(fname))
            total_size += os.path.getsize(fname)
            os.remove(fname)

    if not forked:
        for fname in filepaths:
            if os.path.isdir(fname):
                shutil.rmtree(fname)

    return total_size


def remove_files_parallel(filepaths):
    total_size = 0

    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
        futures = []
        for fname in filepaths:
            if not os.path.exists(fname):
                continue
            future = pool.submit(remove_files, [fname], True)
            futures.append(future)

        for future in concurrent.futures.as_completed(futures):
            total_size += future.result()

    return total_size


def check_seqid(seqid):
    taxid = None
    match = re.match(r"(?:^|\|)kraken:taxid\|(\d+)", seqid)
    if match:
        taxid = match.group(1)
    elif re.match(r"^\d+$", seqid):
        taxid = seqid
    if not taxid:
        match = re.match(r"(?:^|\|)([A-Z]+_?[A-Z0-9]+)(?:\||\b|\.)", seqid)
        if match:
            taxid = match.group(1)
    return taxid


def hash_string(string):
    md5 = hashlib.md5()
    md5.update(string.encode())
    return md5.hexdigest()


def hash_file(filename, buf_size=8192):
    LOG.info("Calculating MD5 sum for {}\n".format(filename))
    md5 = hashlib.md5()
    with open(filename, "rb") as in_file:
        while True:
            data = in_file.read(buf_size)
            if not data:
                break
            md5.update(data)
    digest = md5.hexdigest()
    LOG.info("MD5 sum of {} is {}\n".format(filename, digest))
    return digest


# This function is part of the Kraken 2 taxonomic sequence
# classification system.
#
# Reads multi-FASTA input and examines each sequence header.  Headers are
# OK if a taxonomy ID is found (as either the entire sequence ID or as part
# of a "kraken:taxid" token), or if something looking like an accession
# number is found.  Not "OK" headers will are fatal errors unless "lenient"
# is used.
#
# Each sequence header results in a line with three tab-separated values;
# the first indicating whether third column is the taxonomy ID ("TAXID") or
# an accession number ("ACCNUM") for the sequence ID listed in the second
# column.
#
def scan_fasta_file(
    in_file, out_file, lenient=False, sequence_to_url=None
):
    LOG.info("Generating prelim_map.txt for {}.\n".format(in_file.name))
    iterator = in_file
    iterator_is_dict = False
    if type(sequence_to_url) is dict:
        iterator = sequence_to_url
        iterator_is_dict = True
    for line in iterator:
        if not line.startswith(">"):
            continue
        remote_filepath = sequence_to_url
        if iterator_is_dict:
            remote_filepath = sequence_to_url[line]
        for match in re.finditer(r"(?:^>|\x01)(\S+)(?: (.*))?", line):
            seqid = match.group(1)
            taxid = check_seqid(seqid)
            comment = match.group(2) or ""
            if not taxid:
                if lenient:
                    continue
                else:
                    sys.exit(1)
            if re.match(r"^\d+$", taxid):
                out_file.write(
                    "TAXID\t{:s}\t{:s}\t{:s}\t{:s}\n".format(
                        seqid, taxid, comment, remote_filepath
                    )
                )
            else:
                out_file.write(
                    "ACCNUM\t{:s}\t{:s}\t{:s}\t{:s}\n".format(
                        seqid, taxid, comment, remote_filepath
                    )
                )
    LOG.info(
        "Finished generating prelim_map.txt for {}.\n".format(in_file.name)
    )


# This function is part of the Kraken 2 taxonomic sequence
# classification system.
#
# Looks up accession numbers and reports associated taxonomy IDs
#
# `lookup_list_file` is 1 2-column TSV file w/ sequence IDs and
# accession numbers, and `accession_map_files` is a list of
# accession2taxid files from NCBI.  Output is tab-delimited lines,
# with sequence IDs in first column and taxonomy IDs in second.
#
def lookup_accession_numbers(
    lookup_list_filename, out_filename, *accession_map_files
):
    target_lists = {}
    with open(lookup_list_filename, "r") as f:
        for line in f:
            line = line.strip()
            seqid, acc_num = line.split("\t")
            if acc_num in target_lists:
                target_lists[acc_num].append(seqid)
            else:
                target_lists[acc_num] = [seqid]
    initial_target_count = len(target_lists)
    with open(out_filename, "a") as out_file:
        for filename in accession_map_files:
            with open(filename, "r") as in_file:
                in_file.readline()  # discard header line
                line_count = 0
                for line in in_file:
                    line_count += 1
                    line = line.strip()
                    split = line.split("\t")
                    if len(split) != 4:
                        LOG.warning(
                            "{}:{}-'{}' contains fewer than 4 fields\n"
                            .format(filename, line_count, line)
                        )
                        continue
                    accession, with_version, taxid, gi = split
                    if accession in target_lists:
                        lst = target_lists[accession]
                        del target_lists[accession]
                        for seqid in lst:
                            out_file.write(seqid + "\t" + taxid + "\n")
                        if len(target_lists) == 0:
                            break
            if len(target_lists) == 0:
                break
    if target_lists:
        LOG.warning(
            "{}/{} accession numbers remain unmapped, "
            "see unmapped_accessions.txt in {} directory\n"
            .format(len(target_lists), initial_target_count,
                    os.path.abspath(os.curdir))
        )
        with open("unmapped_accessions.txt", "w") as f:
            for k in target_lists:
                f.write(k + "\n")


def spawn_masking_subprocess(output_file, threads, protein=False):
    masking_binary = "segmasker" if protein else "k2mask"
    if "MASKER" in os.environ:
        masking_binary = os.environ["MASKER"]
    masking_binary = find_kraken2_binary(masking_binary)

    argv = masking_binary + " -outfmt fasta | sed -e '/^>/!s/[a-z]/x/g'"
    if masking_binary.find("k2mask") >= 0:
        # k2mask can run multithreaded
        argv = masking_binary + " -outfmt fasta -threads {} -r x".format(
            threads
        )

    cwd = os.path.dirname(os.path.abspath(output_file.name))
    p = subprocess.Popen(
        argv, shell=True, cwd=cwd,
        stdin=subprocess.PIPE, stdout=output_file
    )

    return p


# Mask low complexity sequences in the database
def mask_files(input_filenames, output_filename, threads, protein=False):
    with open(output_filename, "wb") as fout:
        masker = spawn_masking_subprocess(fout, threads, protein)
        # number_of_files = len(input_filenames)
        for i, input_filename in enumerate(input_filenames):
            library_name = os.path.basename(os.getcwd())
            if library_name == "added":
                LOG.info(
                    "Masking low-complexity regions of added "
                    "library {}\n".format(input_filename)
                )
            else:
                LOG.info(
                    "Masking low-complexity regions of downloaded "
                    "library {:s}\n".format(library_name)
                )
            with open(input_filename, "rb") as fin:
                shutil.copyfileobj(fin, masker.stdin)
                # masker(fin, i + 1 == number_of_files)
        masker.stdin.close()
        if masker.wait() != 0:
            LOG.error("Error while masking {}\n".format(input_filename))


def add_file(args, filename, hashes):
    already_added = False
    filehash = None
    if filename in hashes:
        already_added = True

    filehash = hashes.get(filename) or hash_file(filename)
    destination = os.path.basename(filename)
    ext = ".faa" if args.protein else ".fna"
    base, _ = os.path.splitext(destination)
    destination = base + "_" + filehash + ext
    if already_added:
        LOG.info(
            "Already added " + filename + " to library. "
            "Please remove the entry from added.md5 if this"
            " is not the case.\n"
        )
        return (filename, filehash, destination)

    LOG.info("Adding " + filename + " to library " + args.db + "\n")
    prelim_map_filename = "prelim_map_" + filehash + ".txt"
    with open(prelim_map_filename, mode="a") as out_file:
        with open(filename, "r") as in_file:
            scan_fasta_file(
                in_file, out_file, lenient=True, sequence_to_url=filename
            )
        shutil.copyfile(filename, destination)

    if not args.no_masking:
        mask_files(
            [destination],
            destination + ".masked",
            threads=args.masker_threads,
            protein=args.protein,
        )
        shutil.move(destination + ".masked", destination)
    LOG.info("Added " + filename + " to library " + args.db + "\n")

    return (filename, filehash, destination)


def add_to_library(args):
    if not os.path.isdir(args.db):
        LOG.error("Invalid database: {:s}\n".format(args.db))
        sys.exit(1)
    library_pathname = os.path.join(args.db, "library")
    added_pathname = os.path.join(library_pathname, "added")
    os.makedirs(added_pathname, exist_ok=True)
    args.files = [os.path.abspath(f) for f in args.files]
    os.chdir(added_pathname)
    hashes = {}
    if os.path.exists("added.md5"):
        with open("added.md5", "r") as in_file:
            hashes = dict([line.split()[:2] for line in in_file.readlines()])
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=args.threads
    ) as pool:
        futures = []
        files = map(lambda f: glob.glob(f, recursive=True), args.files)
        for filename in itertools.chain(*files):
            future = pool.submit(add_file, args, filename, hashes)
            if future_raised_exception(future):
                LOG.error(
                    "Error while adding file to library\n"
                )
                raise future.exception()
            futures.append(future)
        with open("added.md5", "a") as out_file:
            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                (filename, filehash, destination) = result
                out_file.write(
                    filename + "\t" + filehash + "\t" + destination + "\n"
                )


def make_manifest_from_assembly_summary(
    assembly_summary_file, is_protein=False
):
    suffix = "_protein.faa.gz" if is_protein else "_genomic.fna.gz"
    manifest_to_taxid = {}
    for line in assembly_summary_file:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        taxid, asm_level, ftp_path = fields[5], fields[11], fields[19]
        if not re.match("Complete Genome|Chromosome", asm_level):
            continue
        if ftp_path == "na":
            continue
        remote_path = ftp_path + "/" + os.path.basename(ftp_path) + suffix
        url_components = urllib.parse.urlsplit(remote_path)
        local_path = url_components.path.replace("/", "", 1)
        manifest_to_taxid[local_path] = taxid
    with open("manifest.txt", "w") as f:
        for k in manifest_to_taxid:
            f.write(k + "\n")
    return manifest_to_taxid


def assign_taxids(args, filepath, manifest_to_taxid,
                  accession_to_taxid={}, filepath_to_url={}):
    sequences_added = 0
    ch_added = 0
    # taxid = manifest_to_taxid[filepath]
    out_filepath = ""
    if filepath.endswith(".gz"):
        out_filepath = os.path.splitext(filepath)[0]
    else:
        out_filepath = filepath + ".tmp"
    masker = None
    sequence_to_url = {}
    with open(out_filepath, "w") as out_file:
        if not args.no_masking:
            masker = spawn_masking_subprocess(
                out_file, args.masker_threads, False
            )
        opener = open
        if filepath.endswith(".gz"):
            opener = gzip.open
        with opener(filepath, "rt") as in_file:
            while True:
                line = in_file.readline()
                if line == "":
                    break
                if line.startswith(">"):
                    taxid = manifest_to_taxid[filepath]
                    if not taxid:
                        match = re.search(r"GC[AF]_[0-9]{9}\.\d+", line)
                        if not match or match.group(0) not in accession_to_taxid:
                            LOG.error(
                                "Unable to assign taxid to sequence {} in file {}\n"
                                .format(line, in_file.name)
                            )
                            sys.exit(1)
                        taxid = accession_to_taxid[match.group(0)]
                    line = line.replace(">", ">kraken:taxid|" + taxid + "|", 1)
                    sequence_to_url[line] = filepath_to_url[filepath]
                    sequences_added += 1
                else:
                    ch_added += len(line) - 1
                if not masker:
                    out_file.write(line)
                else:
                    masker.stdin.write(line.encode())
                taxid = ""
        if out_filepath.endswith(".tmp"):
            shutil.move(out_filepath, filepath)
        if masker:
            masker.stdin.close()
            masker.wait()

    return (sequences_added, ch_added, sequence_to_url)


def download_dataset_by_project(args, endpoint, identifiers):
    library_pathname = os.path.join(args.db, "library")
    os.makedirs(library_pathname, exist_ok=True)
    os.chdir(library_pathname)
    oldwd = os.path.curdir

    for identifier in identifiers:
        dirname = identifier.lower().replace(" ", "_")
        os.makedirs(dirname, exist_ok=True)
        os.chdir(dirname)
        download_and_process_accessions(args, endpoint, identifier)
        os.chdir(oldwd)


def download_and_process_accessions(args, endpoint, identifier):
    api = http.client.HTTPSConnection(NCBI_REST_API)
    identifier = urllib.parse.quote(identifier)
    accession_to_taxid = {}
    builder = NCBI_URI_Builder(endpoint, identifier, "dataset_report")
    builder.set_filters_from_args(args)
    builder = builder.page_size(100)
    old_page_token = ""

    while True:
        api.request("GET", builder.build())
        response = api.getresponse()
        results = response.readlines()[0]
        results = json.loads(results)
        response.close()

        if not results:
            LOG.error(
                "Could not find any accessions matching the query: {}\n"
                .format(identifier)
            )
            sys.exit(1)

        for report in results["reports"]:
            accession_to_taxid[report["accession"]] =\
                report["organism"]["tax_id"]
        if "next_page_token" in results \
           and results["next_page_token"] != old_page_token:
            builder = builder.page_token(results["next_page_token"])
            old_page_token = results["next_page_token"]
        else:
            break

    api.close()
    LOG.info(
        "Found {} accession(s) associated with {}\n"
        .format(len(accession_to_taxid), identifier)
    )

    accessions = list(accession_to_taxid.keys())
    accession_to_url = map_accessions_to_url_parallel(args, accessions)

    filepath_to_taxid = {}
    with open("manifest.txt", "w") as fout:
        for accession, url in accession_to_url.items():
            # exclude the leading /
            filepath = urllib.parse.urlparse(url).path[1:]
            fout.write(filepath + "\n")
            taxid = accession_to_taxid[accession]
            filepath_to_taxid[filepath] = str(taxid)

    download_files_from_manifest(NCBI_SERVER, args.threads)
    # filepath_to_taxid = download_accessions(args, accession_to_taxid)
    filepath_to_url = {}
    for filepath in filepath_to_taxid.keys():
        accession = re.search(r"GC[AF]_\d{9}\.\d+", filepath).group()
        filepath_to_url[filepath] = accession_to_url[accession]
    sequence_to_url = assign_taxid_to_sequences(
        args, filepath_to_taxid, filepath_to_url=filepath_to_url
    )

    library_filename = "library.faa" if args.protein else "library.fna"
    with open(library_filename, "r") as in_file:
        with open("prelim_map.txt", "w") as out_file:
            out_file.write("# prelim_map for " + args.library + "\n")
            scan_fasta_file(
                in_file,
                out_file,
                sequence_to_url=sequence_to_url,
            )


def download_accessions(args, accession_to_taxid):
    accessions = sorted(list(accession_to_taxid.keys()))
    filepath_to_taxid = {}
    zip_file_list = []

    zip_filename_to_accessions = {}
    extension = ".faa" if args.protein else ".fna"
    if os.path.exists("zips"):
        downloaded_accessions = {}
        for zip_filename in os.listdir("zips"):
            zip_filename = os.path.join("zips", zip_filename)
            if zip_filename.endswith(".zip"):
                with zipfile.ZipFile(zip_filename) as zip:
                    for filename in zip.namelist():
                        match = re.search(r"GC[AF]_\d{9}\.\d+", filename)
                        if match and filename.endswith(extension):
                            accession = match.group()
                            downloaded_accessions[accession] = zip_filename
                zip_file_list.append(zip_filename)

        if downloaded_accessions:
            unfetched_accessions = []
            for accession in accessions:
                if accession in downloaded_accessions:
                    LOG.info(
                        "Already downloaded accession: {}, skipping\n"
                        .format(accession)
                    )
                    zip_filename = downloaded_accessions[accession]
                    if zip_filename in zip_filename_to_accessions:
                        zip_filename_to_accessions[zip_filename].append(accession)
                    else:
                        zip_filename_to_accessions[zip_filename] = [accession]
                else:
                    unfetched_accessions.append(accession)
            accessions = unfetched_accessions

    partitions = []
    if accessions:
        number_of_partitions = math.ceil(len(accessions) / 400)
        number_of_partitions = max(args.threads, number_of_partitions)
        partitions = partition_list(accessions, number_of_partitions)

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as pool:
        download_futures = []
        for partition in partitions:
            download_futures.append(
                pool.submit(download_zip_from_ncbi, args, partition)
            )
        unzip_futures = []
        for future in concurrent.futures.as_completed(download_futures):
            unzip_futures.append(
                pool.submit(
                    extract_fastas_from_zip_file,
                    future.result(), args.protein
                )
            )
        for zip_filename, accessions in zip_filename_to_accessions.items():
            unzip_futures.append(
                pool.submit(
                    extract_fastas_from_zip_file,
                    zip_filename, args.protein,
                    accessions
                )
            )
        for future in concurrent.futures.as_completed(unzip_futures):
            accession_to_filepath = future.result()
            for accession, filepath in accession_to_filepath.items():
                filepath_to_taxid[filepath] = str(accession_to_taxid[accession])

    shutil.rmtree("ncbi_dataset", ignore_errors=True)

    return filepath_to_taxid


def extract_fastas_from_zip_file(filename, protein, accessions_to_extract=None):
    accession_to_filepath = {}

    suffix = "_protein.faa" if protein else "_genomic.fna"
    with zipfile.ZipFile(filename) as zip:
        for entry in zip.namelist():
            if entry.endswith(suffix):
                dir_components = ["genomes", "all"]
                filename = os.path.basename(entry)
                dirname = filename.replace(suffix, "")
                accession = re.search(r"GC[AF]_\d{9}\.\d+", filename).group()
                if accessions_to_extract and\
                   accession not in accessions_to_extract:
                    continue
                modified_accession = accession.split(".")[0].replace("_", "")
                dir_components.extend(partition_list(modified_accession, 4))
                dir_components.append(dirname)
                dir_components.append(filename)

                filepath = os.path.join("", *dir_components)
                entry_size = zip.getinfo(entry).file_size
                if os.path.exists(filepath) and os.stat(filepath).st_size == entry_size:
                    LOG.info(
                        "Already extracted {} from {}\n"
                        .format(entry, os.path.abspath(filename))
                    )
                else:
                    LOG.info(
                        "Extracting {} from {}\n".format(entry, os.path.abspath(filename))
                    )
                    zip.extract(entry, os.path.curdir)
                    os.makedirs(os.path.dirname(filepath), exist_ok=True)
                    LOG.debug(
                        "Moving {} to {}\n".format(
                            os.path.abspath(entry), os.path.abspath(filepath)
                        )
                    )
                    shutil.move(entry, filepath)
                accession_to_filepath[accession] = filepath

    return accession_to_filepath


def download_zip_from_ncbi(args, accessions):
    api = http.client.HTTPSConnection(NCBI_REST_API)
    accessions = ",".join(accessions)
    md5 = hash_string(accessions)
    os.makedirs("zips", exist_ok=True)
    filename = os.path.join("zips", md5 + ".zip")
    tmp_filename = filename + ".tmp"
    if os.path.exists(os.path.join("zip", filename)):
        LOG.info(
            "Already downloaded {} from NCBI which contains accessions: {}\n"
            .format(os.path.basename(filename), accessions)
        )
        return filename
    LOG.info(
        "Downloading {} from NCBI containing the following accessions: {}\n"
        .format(os.path.basename(filename), accessions)
    )
    accessions = urllib.parse.quote(accessions)
    annotation_type = "PROT_FASTA" if args.protein else "GENOME_FASTA"
    builder = NCBI_URI_Builder("accession", accessions, "download")
    builder = builder.include_annotation_type(annotation_type)

    api.request("GET", builder.build())
    res = api.getresponse()

    with open(tmp_filename, "wb") as fout:
        shutil.copyfileobj(res, fout)

    res.close()
    shutil.move(tmp_filename, filename)
    LOG.info("Saved {} to {}\n".format(
        os.path.basename(filename),
        os.path.abspath(filename)
    ))

    return filename


def partition_list(list, num_partitions):
    partitions = []
    step = math.ceil(len(list) / num_partitions)
    for i in range(0, len(list), step):
        partitions.append(list[i:i+step])

    return partitions


def map_accessions_to_url(accessions, protein=False):
    api = http.client.HTTPSConnection(NCBI_REST_API)
    conn = http.client.HTTPSConnection(NCBI_SERVER)

    size = len(accessions)
    start = 0
    step = 20
    stop = min(size, step)
    base_uri = "/datasets/v2alpha/genome/accession/{}/links"
    accession_to_url = {}

    while True:
        query = ",".join(accessions[start:stop])
        query = urllib.parse.quote(query)
        api.request("GET", base_uri.format(query))
        req = api.getresponse()
        results = req.readlines()[0]
        req.close()
        results = json.loads(results)

        for entry in results["assembly_links"]:
            if entry["assembly_link_type"] == "FTP_LINK":
                filepath = get_download_path(
                    conn, entry["resource_link"], protein
                )
                accession_to_url[entry["accession"]] = url_join(NCBI_SERVER, path=filepath)
        start, stop = stop, min(size, stop + step)

        if start == stop:
            break
    api.close()

    return accession_to_url


def map_accessions_to_url_parallel(args, accessions):
    futures = []
    accession_to_url = {}
    LOG.info("Fetching download links for accessions\n")
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as pool:
        for partition in partition_list(accessions, args.threads):
            future = pool.submit(
                map_accessions_to_url, partition, args.protein
            )
            futures.append(future)
    for future in concurrent.futures.as_completed(futures):
        accession_to_url.update(future.result())
    LOG.info("Finished fetching download links for accessions\n")

    return accession_to_url


def get_download_path(conn, resource_link, protein):
    dir = os.path.basename(resource_link)
    filename = dir + "_genomic.fna.gz"
    resource_link += "/"
    if protein:
        filename = dir + "_protein.faa.gz"
    resource_link = urllib.parse.urljoin(resource_link, filename)
    path = urllib.parse.urlparse(resource_link).path
    # conn.request("HEAD", path)
    # response = conn.getresponse()
    # if response.status != 200:
    #     LOG.error("Unable to find ...")
    #     sys.exit(1)
    # response.close()

    return path


def assign_taxid_to_sequences(args,
                              manifest_to_taxid,
                              accession_to_taxid={},
                              filepath_to_url={}):
    if args.no_masking:
        LOG.info("Assigning taxonomic IDs to sequences\n")
    else:
        LOG.info(
            "Assigning taxonomic IDs and masking sequences\n"
        )
    library_filename = "library.faa" if args.protein else "library.fna"
    sequence_to_url = {}
    projects_added = 0
    total_projects = len(manifest_to_taxid)
    sequences_added = 0
    ch_added = 0
    ch = "aa" if args.protein else "bp"

    out_line = progress_line(
        projects_added, total_projects, sequences_added, ch_added, ch
    )

    LOG.debug("{:s}\r".format(out_line))
    filepaths = sorted(manifest_to_taxid)
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as pool:
        futures = []
        max_out_line_len = 0
        for filepath in filepaths:
            future = pool.submit(
                assign_taxids, args, filepath,
                manifest_to_taxid, accession_to_taxid,
                filepath_to_url
            )
            if future_raised_exception(future):
                LOG.error(
                    "Error encountered while assigning tax IDs\n"
                )
                raise future.exception()
            futures.append(future)
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            sequences_added += result[0]
            ch_added += result[1]
            projects_added += 1
            sequence_to_url.update(result[2])
            out_line = progress_line(
                projects_added, total_projects, sequences_added, ch_added, ch
            )
            max_out_line_len = max(len(out_line), max_out_line_len)
            padding = " " * (max_out_line_len - len(out_line))
            LOG.debug("{:s}\r".format(out_line + padding))

    if args.no_masking:
        LOG.info("Finished assigning taxonomic IDs to sequences\n")
    else:
        LOG.info("Finished assigning taxonomic IDs and masking sequences\n")

    LOG.info("Generating {:s}\n".format(library_filename))
    with open(library_filename, "w") as out_file:
        for filepath in filepaths:
            if filepath.endswith(".gz"):
                filepath = os.path.splitext(filepath)[0]
            with open(filepath, "r") as in_file:
                shutil.copyfileobj(in_file, out_file)

    LOG.info("Finished generating {:s}\n".format(library_filename))

    return sequence_to_url


def progress_line(projects, total_projects, seqs, chars, ch):
    line = "Processed "
    if projects == total_projects:
        line += str(projects)
    else:
        line += "{:d}/{:d}".format(projects, total_projects)
    line += " project(s), {:d} sequence(s), ".format(seqs)
    prefix = None
    for p in ["k", "M", "G", "T", "P", "E"]:
        if chars >= 1024:
            prefix = p
            chars /= 1024
        else:
            break
    if prefix:
        line += "{:.2f} {:s}{:s}".format(chars, prefix, ch)
    else:
        line += "{:.2f} {:s}".format(chars, ch)
    return line


# The following three functions have dummy return values. This is
# so that we can check whether a future stopped running as a
# result of an exception, so that we can stop the main process
# early.
# def decompress_files(compressed_filenames, out_filename=None, buf_size=8192):
#     if isinstance(compressed_filenames, str):
#         compressed_filenames = [compressed_filenames]
#     if out_filename:
#         if os.path.exists(out_filename + ".tmp"):
#             os.remove(out_filename + ".tmp")
#         with open(out_filename + ".tmp", "ab") as out_file:
#             for filename in compressed_filenames:
#                 with gzip.open(filename) as gz:
#                     decompress_file(gz, out_file)
#             os.rename(out_filename + ".tmp", out_filename)
#     else:
#         for filename in compressed_filenames:
#             out_filename, ext = os.path.splitext(filename)
#             if os.path.exists(out_filename + ".tmp"):
#                 os.remove(out_filename + ".tmp")
#             with gzip.open(filename) as gz:
#                 with open(out_filename + ".tmp", "wb") as out:
#                     decompress_file(gz, out, buf_size)
#             os.rename(out_filename + ".tmp", out_filename)

#     return True


def decompress_files(compressed_filenames, out_filename=None, buf_size=8192):
    if isinstance(compressed_filenames, str):
        compressed_filenames = [compressed_filenames]
    if out_filename:
        if os.path.exists(out_filename + ".tmp"):
            os.remove(out_filename + ".tmp")
        with open(out_filename + ".tmp", "ab") as out_file:
            for filename in compressed_filenames:
                with open(filename, "rb") as gz:
                    decompress_file(gz, out_file)
            os.rename(out_filename + ".tmp", out_filename)
    else:
        for filename in compressed_filenames:
            out_filename, ext = os.path.splitext(filename)
            if os.path.exists(out_filename + ".tmp"):
                os.remove(out_filename + ".tmp")
            with open(filename, "rb") as gz:
                with open(out_filename + ".tmp", "wb") as out:
                    decompress_file(gz, out, buf_size)
            os.rename(out_filename + ".tmp", out_filename)

    return True


def decompress_file(in_file, out_file, buf_size=8129):
    LOG.info(
        "Decompressing {:s}\n".format(os.path.join(os.getcwd(), in_file.name))
    )
    inflator = zlib.decompressobj(15 + 32)
    while True:
        data = in_file.read(buf_size)
        if not data:
            break
        out_file.write(inflator.decompress(data))
    # shutil.copyfileobj(in_file, out_file, buf_size)
    LOG.info(
        "Finished decompressing {:s}\n".format(
            os.path.join(os.getcwd(), in_file.name)
        )
    )

    return True


def decompress_and_mask(filepath, masker_threads):
    out_filepath = os.path.splitext(filepath)[0]
    with open(out_filepath, "w") as out_file:
        masker = spawn_masking_subprocess(out_file, masker_threads)
        with gzip.open(filepath, "rb") as in_file:
            decompress_file(in_file, masker.stdin)
        masker.stdin.close()
        masker.wait()

    return True


def download_log(filename, total_size=None):
    pb = None
    current_size = 0

    def inner(block_number, read_size, size):
        nonlocal pb, current_size, total_size
        if not pb:
            pb = ProgressBar(total_size or size)
        current_size += read_size
        pb.progress(current_size)
        LOG.debug(
            "{:s} {: >10s}\r".format(pb.get_bar(), format_bytes(current_size))
        )

    return inner


def http_download_file(url, local_name=None, call_back=None):
    if not local_name:
        local_name = urllib.parse.urlparse(url).path.split("/")[-1]
    else:
        local_name = os.path.abspath(local_name)
        os.makedirs(os.path.dirname(local_name), exist_ok=True)
    with urllib.request.urlopen(url) as conn:
        remote_size = int(conn.headers["Content-Length"])
        local_size = (
            os.stat(local_name).st_size if os.path.exists(local_name) else 0
        )
        if local_size == remote_size:
            LOG.info(
                "Already downloaded {:s}\n".format(get_abs_path(local_name))
            )
            return

    LOG.info("Beginning download of {:s}\n".format(url))
    urllib.request.urlretrieve(
        url, local_name, reporthook=(call_back or download_log(local_name))
    )
    clear_console_line()
    LOG.info("Saved {:s} to {:s}\n".format(local_name, os.getcwd()))


def http_download_file2(server, urls, save_to=None, md5sums=None):
    files_downloaded = 0
    conn = None
    if isinstance(server, str):
        conn = http.client.HTTPSConnection(server, timeout=60)
    else:
        conn = server
        server = conn.host
    md5 = md5sums if md5sums else {}
    for url in urls:
        url = url.strip()
        filename = os.path.basename(url)
        local_name = os.path.abspath(url)
        if save_to:
            local_name = os.path.join(save_to, filename)
        local_directory = os.path.dirname(local_name)
        os.makedirs(os.path.dirname(local_name), exist_ok=True)

        try:
            if os.path.exists(local_name):
                if filename not in md5:
                    checksums = []
                    remote_dirname = os.path.dirname(url)
                    for md5_filename in ["md5checksums.txt", filename + ".md5", "MD5SUM.txt"]:
                        # Check if file exists before trying to download.
                        # This avoids NCBI sending BadStatusLine when
                        # making the request.
                        conn.request(
                            "HEAD", "/" + remote_dirname + "/" + md5_filename
                        )
                        response = conn.getresponse()
                        if response.status == 200:
                            response.close()
                            # If we have found the file, then go ahead
                            # and download.
                            conn.request(
                                "GET", "/" + remote_dirname + "/" + md5_filename
                            )
                            response = conn.getresponse()
                            checksums = response.readlines()
                            response.close()
                            break
                        else:
                            response.close()
                    if len(checksums) > 0:
                        for checksum in checksums:
                            (md5sum, remote_filename) = checksum.split()
                            remote_filename = os.path.basename(
                                remote_filename.decode()
                            )
                            md5[remote_filename] = os.path.basename(md5sum.decode())
                if filename in md5 and md5[filename] == hash_file(local_name):
                    LOG.info(
                        "Already downloaded {:s}\n".format(
                            urllib.parse.urljoin(server, url)
                        )
                    )
                    # Server can potentially end the connection while we
                    # waiting for the md5 hash to be computed. We have
                    # to reconnect to avoid failures when trying to
                    # retrieve the file.
                    conn.connect()
                    continue
            LOG.info("Beginning download of {:s}\n".format(server + "/" + url))
            with open(local_name, "wb") as out_file:
                conn.request("GET", "/" + url)
                response = conn.getresponse()
                if response.status == 200:
                    shutil.copyfileobj(response, out_file, 8192)
                elif response.status == 404:
                    LOG.warning(
                        "Cannot find file: {}.\n"
                        "Please report this issue to NCBI.\n".format(url)
                    )
                else:
                    LOG.error(
                        "Error downloading file: {}.\n"
                        "Reason: {}\n".format(url, response.reason)
                    )
                    response.read()
                    sys.exit(1)
                response.close()

        except http.client.RemoteDisconnected:
            conn.connect()
            LOG.warning("Unable to download " + url + ", will try again\n")
            urls.append(url)
            continue

        LOG.info("Saved {:s} to {:s}\n".format(filename, local_directory))
        files_downloaded += 1

    return files_downloaded


def make_file_filter(file_handle, regex):
    def inner(listing):
        path = listing.split()[-1]
        if path.endswith(regex):
            file_handle.write(path + "\n")

    return inner


def move(src, dst):
    src = os.path.abspath(src)
    dst = os.path.abspath(dst)
    if os.path.isfile(src) and os.path.isdir(dst):
        dst = os.path.join(dst, os.path.basename(src))
    shutil.move(src, dst)


def get_manifest_and_md5sums(server, remote_directory, regex):
    ftp = ftplib.FTP(server)
    ftp.login()

    sstream = io.StringIO()
    ftp.cwd(remote_directory)
    ftp.retrlines("LIST", callback=make_file_filter(sstream, regex))

    with open("manifest.txt", "w") as out:
        for line in sstream.getvalue().split():
            out.write(urllib.parse.urljoin(remote_directory, line) + "\n")

    sstream.truncate(0)
    sstream.seek(0)
    ftp.cwd("/refseq/release/release-catalog")
    ftp.retrlines("LIST", callback=make_file_filter(sstream, "installed"))
    install_file = sstream.getvalue().strip()

    bstream = io.BytesIO()

    ftp.retrbinary("RETR " + install_file, callback=bstream.write)
    ftp.close()

    md5sums = {}
    for line in bstream.getvalue().split(b"\n"):
        if line.find(b"plasmid") == -1:
            continue
        (md5, filename) = line.split()
        md5sums[filename.decode()] = md5.decode()
    return md5sums


def download_files_from_manifest(
    server,
    threads=1,
    manifest_filename="manifest.txt",
    filepath_to_taxid_table=None,
    md5sums=None,
):
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as pool:
        with open(manifest_filename, "r") as f:
            filepaths = f.readlines()
            # We try to reduce the risk of a single thread downloading
            # many large files.
            random.shuffle(filepaths)
            partitions = []
            futures = []
            step = len(filepaths) // threads
            if step == 0:
                step = 1
            for i in range(0, len(filepaths), step):
                partitions.append(filepaths[i:i+step])
            for partition in partitions:
                future = pool.submit(
                    http_download_file2, server, partition, None, md5sums
                )
                if future_raised_exception(future):
                    LOG.error(
                        "Error encounted while trying to download files\n"
                    )
                    raise future.exception()
                futures.append(future)
            concurrent.futures.wait(futures)


def download_taxonomy(args):
    taxonomy_path = os.path.join(args.db, "taxonomy")
    os.makedirs(taxonomy_path, exist_ok=True)
    os.chdir(taxonomy_path)
    futures = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:
        maps = set()
        if not args.skip_maps:
            if not args.protein:
                for subsection in ["gb", "wgs"]:
                    filename = "pub/taxonomy/accession2taxid/"
                    filename += "nucl_" + subsection + ".accession2taxid.gz"
                    maps.add(os.path.basename(filename))
                    future = pool.submit(
                        http_download_file2,
                        NCBI_SERVER,
                        [filename],
                        save_to=os.path.abspath(os.curdir),
                    )
                    if future_raised_exception(future):
                        LOG.error(
                            "Error encountered while downloading file\n"
                        )
                        raise future.exception()
                    futures.append(future)
            else:
                filename = "/pub/taxonomy/accession2taxid/"
                filename += "prot.accession2taxid.gz"
                future = pool.submit(
                    http_download_file2,
                    NCBI_SERVER,
                    [filename],
                    save_to=os.path.abspath(os.curdir),
                )
                if future_raised_exception(future):
                    LOG.error(
                        "Error encountered while downloading file\n"
                    )
                    raise future.exception()
                futures.append(future)
                maps.add(os.path.basename(filename))

            for future in concurrent.futures.as_completed(futures):
                download_count = future.result()
                maps_downloaded = set(glob.glob("*.accession2taxid.gz"))
                map_to_decompress = maps_downloaded.intersection(maps)
                for gz_filename in map_to_decompress:
                    if not os.path.exists(gz_filename[:-3])\
                       or download_count != 0:
                        future = pool.submit(decompress_files, [gz_filename])
                        if future_raised_exception(future):
                            LOG.error(
                                "Error encountered while decompressing file\n"
                            )
                            raise future.exception()
                        futures.append(future)
                    maps.remove(gz_filename)

    LOG.info("Downloading taxonomy tree data\n")
    filename = "pub/taxonomy/taxdump.tar.gz"
    http_download_file2(
        NCBI_SERVER, [filename], save_to=os.path.abspath(os.curdir)
    )
    LOG.info("Untarring taxonomy tree data\n")
    with tarfile.open("taxdump.tar.gz", "r:gz") as tar:
        tar.extractall()
    LOG.info("Finished Untarring taxonomy tree data")


def download_gtdb_taxonomy(args, files, md5s):
    taxonomy_path = os.path.join(args.db, "taxonomy")
    os.makedirs(taxonomy_path, exist_ok=True)
    os.chdir(taxonomy_path)

    LOG.info("Dowloading GTDB taxonomy for bacteria and archaea\n")
    http_download_file2(
        GTDB_SERVER,
        files,
        save_to=os.path.abspath(os.curdir), md5sums=md5s
    )
    LOG.info("Finished downloading GTDB taxonomy for bacteria and archaea\n")


def build_gtdb_taxonomy(in_file):
    rank_codes = {
        "d": "domain",
        "p": "phylum",
        "c": "class",
        "o": "order",
        "f": "family",
        "g": "genus",
        "s": "species",
    }
    accession_map = {}
    seen_it = collections.defaultdict(int)
    child_data = collections.defaultdict(lambda: collections.defaultdict(int))
    for line in in_file:
        line = line.strip()
        accession, taxonomy_string = line.split("\t")
        start = accession.find("GCA")
        if start < 0:
            start = accession.find("GCF")
        accession = accession[start:]
        taxonomy_string = re.sub("(;[a-z]__)+$", "", taxonomy_string)
        accession_map[accession] = taxonomy_string
        seen_it[taxonomy_string] += 1
        if seen_it[taxonomy_string] > 1:
            continue
        while True:
            match = re.search("(;[a-z]__[^;]+$)", taxonomy_string)
            if not match:
                break
            level = match.group(1)
            taxonomy_string = re.sub("(;[a-z]__[^;]+$)", "", taxonomy_string)
            key = taxonomy_string + level
            child_data[taxonomy_string][key] += 1
            seen_it[taxonomy_string] += 1
            if seen_it[taxonomy_string] > 1:
                break
        if seen_it[taxonomy_string] == 1:
            child_data["root"][taxonomy_string] += 1

    id_map = {}
    next_node_id = 1
    LOG.info("Generating nodes.dmp and names.dmp\n")
    with open("names.dmp", "w") as names_file:
        with open("nodes.dmp", "w") as nodes_file:
            bfs_queue = [["root", 1]]
            while len(bfs_queue) > 0:
                node, parent_id = bfs_queue.pop()
                display_name = node
                rank = None
                match = re.search("([a-z])__([^;]+)$", node)
                if match:
                    rank = rank_codes[match.group(1)]
                    display_name = match.group(2)
                rank = rank or "no rank"
                node_id, next_node_id = next_node_id, next_node_id + 1
                id_map[node] = node_id
                names_file.write(
                    "{:d}\t|\t{:s}\t|\t-\t|\tscientific name\t|\n".format(
                        node_id, display_name
                    )
                )
                nodes_file.write(
                    "{:d}\t|\t{:d}\t|\t{:s}\t|\t-\t|\n".format(
                        node_id, parent_id, rank
                    )
                )
                children = (
                    sorted([key for key in child_data[node]])
                    if node in child_data
                    else []
                )
                for node in children:
                    bfs_queue.insert(0, [node, node_id])
    with open("gtdb.accession2taxid", "w") as f:
        for accession in sorted([key for key in accession_map]):
            taxid = id_map[accession_map[accession]]
            accession_without_revision = accession.split(".")[0]
            f.write("{:s}\t{:s}\t{:d}\t-\n".format(
                accession_without_revision,
                accession, taxid
            ))


def download_gtdb_genomes(args, remote_filepath, md5s):
    for directory in ["taxonomy", "library"]:
        os.makedirs(directory, exist_ok=True)
    filename = os.path.basename(remote_filepath)
    filepath = os.path.abspath(filename)
    http_download_file2(
        GTDB_SERVER, [remote_filepath],
        save_to=os.curdir, md5sums=md5s
    )
    os.chdir("library")
    accession_to_filepath = {}
    filepaths_without_accession = []
    ext = ".faa" if args.protein else ".fna"
    library = ""
    if filename.endswith(".tar.gz"):
        with tarfile.open(filepath, "r:gz") as tar:
            library = tar.getnames()[0]
            for member in tar.getmembers():
                if member.isfile() and re.search(ext, member.name):
                    if re.search("GC[AF]", member.name):
                        filename = os.path.basename(member.name)
                        accession = re.search(
                            r"(GC[AF]_\d{9}\.\d+)", filename
                        ).group(1)
                        if os.path.exists(member.name)\
                           and member.size == os.stat(member.name).st_size:
                            LOG.info(
                                "Already extracted {}...skipping\n"
                                .format(member.name)
                            )
                        else:
                            LOG.info("Extracting {}\n".format(member.name))
                            tar.extract(member)
                            LOG.info(
                                "Finished extracting {}\n"
                                .format(member.name)
                            )
                        accession_to_filepath[accession] =\
                            os.path.abspath(member.name)
                    else:
                        # We do not check if a file has already been extracted here
                        LOG.info("Extracting {}\n".format(member.name))
                        tar.extract(member)
                        LOG.info(
                            "Finished extracting {}\n".format(member.name)
                        )
                        filepaths_without_accession.append(
                            os.path.abspath(member.name)
                        )
    else:
        library = filename.split(".")[0]
        filepaths_without_accession.append(os.path.abspath(filepath))

    return (library, accession_to_filepath, filepaths_without_accession)


def build_gtdb_database(args):
    db_pathname = os.path.abspath(args.db)
    os.makedirs(db_pathname, exist_ok=True)
    os.chdir(db_pathname)
    md5_url = url_join(GTDB_SERVER, path="/releases/latest/MD5SUM.txt")
    http_download_file(md5_url)
    files_needed = collections.defaultdict(list)
    md5s = {}

    with open("MD5SUM.txt", "r") as in_file:
        # filename_regex = r"genomes|genes|fna"
        filename_regex = "|".join(args.gtdb_files)
        pattern = re.compile(filename_regex)
        candidate_files = []
        if args.protein:
            filename_regex = r"protein|faa"
        for line in in_file:
            md5sum, filepath = line.split()
            filepath = urllib.parse.urljoin("releases/latest/", filepath)
            # remove the release tag since they do not appear in the "latest"
            # file listings
            filepath = re.sub(r"_r\d+", "", filepath)
            if filepath.find("genomic_files") != -1:
                candidate_files.append(os.path.basename(filepath))
            if re.search(r"taxonomy.*\.tsv$", filepath):
                files_needed["taxonomy"].append(filepath)
            elif re.search(pattern, filepath):
                files_needed["fasta"].append(filepath)
            filepath = os.path.basename(filepath)
            md5s[filepath] = md5sum

    if len(files_needed["fasta"]) != len(args.gtdb_files):
        LOG.error("At least one of the files did not match: {}\n"
                  .format(", ".join(args.gtdb_files)))
        LOG.error("Here are a list of candidates:\n{}\n"
                  .format("\n".join(candidate_files)))
        sys.exit(1)
    download_gtdb_taxonomy(args, files_needed["taxonomy"], md5s)
    os.chdir(os.path.join(db_pathname, "taxonomy"))
    LOG.info("Merging Archaea and Bacteria taxonomies\n")
    with open("merged_taxonomy.tsv", "w") as file_out:
        for tax_filename in files_needed["taxonomy"]:
            tax_filename = os.path.basename(tax_filename)
            with open(tax_filename, "r") as file_in:
                shutil.copyfileobj(file_in, file_out)
    LOG.info("Finished merging Archaea and Bacteria taxonomies\n")

    with open("merged_taxonomy.tsv", "r") as in_file:
        build_gtdb_taxonomy(in_file)

    workers = len(files_needed["fasta"])
    futures = []
    accession_to_filepath = {}
    # These files are not tied to a single accession, but
    # instead each sequence in the FASTA has its own accession.

    os.chdir(os.path.join(db_pathname, "taxonomy"))
    accession_to_taxid = {}
    with open("gtdb.accession2taxid", "r") as in_file:
        for line in in_file:
            base_accession, accession, taxid, gi = line.split("\t")
            accession_to_taxid[accession] = taxid

    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
        for remote_filepath in files_needed["fasta"]:
            os.chdir(db_pathname)
            futures.append(pool.submit(
                download_gtdb_genomes, args, remote_filepath, md5s
            ))
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            library = result[0]
            filepaths_without_accessions = []
            filepath_to_url = {}

            for accession, filepath in result[1].items():
                accession_to_filepath[accession] = filepath
                filepath_to_url[filepath] =\
                    url_join(GTDB_SERVER, path=remote_filepath)
            for filepath in result[2]:
                filepaths_without_accessions.append(filepath)
                filepath_to_url[filepath] =\
                    url_join(GTDB_SERVER, path=remote_filepath)

            filepath_to_taxid_table = {}
            library_pathname = os.path.join(db_pathname,
                                            os.path.join("library", library))
            os.makedirs(library_pathname, exist_ok=True)
            os.chdir(library_pathname)
            for accession, filepath in accession_to_filepath.items():
                filepath_to_taxid_table[filepath] =\
                    accession_to_taxid[accession]
            for filepath in filepaths_without_accessions:
                filepath_to_taxid_table[filepath] = ""

            sequence_to_url = assign_taxid_to_sequences(
                args, filepath_to_taxid_table,
                accession_to_taxid, filepath_to_url
            )
            with open("library.fna", "r") as in_file:
                with open("prelim_map.txt", "w") as out_file:
                    out_file.write("# prelim_map for {:s}\n".format(library))
                    scan_fasta_file(
                        in_file,
                        out_file,
                        sequence_to_url=sequence_to_url,
                    )

    os.chdir(db_pathname)
    build_kraken2_db(args)


def download_genomic_library(args):
    library_filename = "library.faa" if args.protein else "library.fna"
    library_pathname = os.path.join(args.db, "library")
    LOG.info("Adding {:s} to {:s}\n".format(args.library, args.db))
    if args.library in [
            "archaea",
            "bacteria",
            "viral",
            "fungi",
            "invertebrate",
            "plant",
            "human",
            "protozoa",
            "vertebrate_mammalian",
            "vertebrate_other"
    ]:
        library_pathname = os.path.join(library_pathname, args.library)
        os.makedirs(library_pathname, exist_ok=True)
        os.chdir(library_pathname)
        try:
            os.remove("assembly_summary.txt")
        except FileNotFoundError:
            pass
        remote_dir_name = args.library
        if args.library == "human":
            remote_dir_name = "vertebrate_mammalian/Homo_sapiens"
        try:
            url = "genomes/refseq/{:s}/assembly_summary.txt".format(
                remote_dir_name
            )
            http_download_file2(
                NCBI_SERVER, [url], save_to=os.path.abspath(os.curdir)
            )
        except urllib.error.URLError:
            LOG.error(
                "Error downloading assembly summary file for {:s}, "
                "exiting\n".format(args.library)
            )
            sys.exit(1)
        if args.library == "human":
            with open("assembly_summary.txt", "r") as f1:
                with open("grc.txt", "w") as f2:
                    for line in f1:
                        if line.find("Genome Reference Consortium"):
                            f2.write(line)
            os.rename("grc.txt", "assembly_summary.txt")
        with open("assembly_summary.txt", "r") as f:
            filepath_to_url = {}
            filepath_to_taxid_table = make_manifest_from_assembly_summary(
                f, args.protein
            )
            for filepath in filepath_to_taxid_table:
                filepath_to_url[filepath] = url_join(NCBI_SERVER, path=filepath)
            download_files_from_manifest(
                NCBI_SERVER,
                args.threads,
                filepath_to_taxid_table=filepath_to_taxid_table,
            )
            sequence_to_url = assign_taxid_to_sequences(
                args, filepath_to_taxid_table,
                filepath_to_url=filepath_to_url
            )
        with open(library_filename, "r") as in_file:
            with open("prelim_map.txt", "w") as out_file:
                out_file.write("# prelim_map for " + args.library + "\n")
                scan_fasta_file(
                    in_file,
                    out_file,
                    sequence_to_url=sequence_to_url,
                )
    elif args.library in ["plasmid", "plastid", "mitochondrion"]:
        library_pathname = os.path.join(args.db, "library")
        library_pathname = os.path.join(library_pathname, args.library)
        library_filename = "library.faa" if args.protein else "library.fna"
        library_filename = os.path.join(library_pathname, library_filename)
        os.makedirs(library_pathname, exist_ok=True)
        os.chdir(library_pathname)
        pat = ".faa.gz" if args.protein else ".fna.gz"
        md5 = get_manifest_and_md5sums(
            NCBI_SERVER, "genomes/refseq/{}/".format(args.library), pat
        )
        download_files_from_manifest(NCBI_SERVER, args.threads, md5sums=md5)
        sequence_to_url = {}
        LOG.info("Generating {}\n".format(library_filename))
        filenames = []
        with open("manifest.txt", "r") as manifest:
            filenames = manifest.readlines()
        filenames = sorted(filenames)
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as pool:
            futures = []
            for filename in filenames:
                filename = filename.strip()
                if not args.no_masking:
                    future = pool.submit(
                        decompress_and_mask, filename, args.masker_threads
                    )
                    if future_raised_exception(future):
                        LOG.error(
                            "Error encountered while decompressing"
                            " or masking files\n"
                        )
                        raise future.exception()
                    futures.append(future)
                else:
                    future = pool.submit(decompress_files, [filename])
                    if future_raised_exception(future):
                        LOG.error(
                            "Error encountered while decompressing files\n"
                        )
                        raise future.exception()
                    futures.append(future)
            result = concurrent.futures.wait(
                futures,
                return_when=concurrent.futures.FIRST_EXCEPTION
            )
            if len(result.not_done) > 0:
                LOG.error(
                    "Encountered error while downloading Plasmid library\n"
                )
                sys.exit(1)
        LOG.info("Generating {}\n".format(library_filename))
        with open(library_filename, "w") as out_file:
            for filename in filenames:
                in_filename = os.path.splitext(filename)[0]
                with open(in_filename, "r") as in_file:
                    for line in in_file:
                        if line.startswith(">"):
                            sequence_to_url[line.strip()] = filename
                        out_file.write(line)

        LOG.info("Finished generating {}\n".format(library_filename))
        with open(library_filename, "r") as in_file:
            with open("prelim_map.txt", "w") as out_file:
                out_file.write("# prelim_map for " + args.library + "\n")
                scan_fasta_file(
                    in_file,
                    out_file,
                    sequence_to_url=sequence_to_url,
                )
    elif args.library in ["nr", "nt"]:
        protein_lib = args.library == "nr"
        if protein_lib and not args.protein:
            LOG.error(
                "{:s} is a protein database, and the Kraken2 database"
                " specified is nucleotide".format(args.library)
            )
            sys.exit(1)
        library_pathname = os.path.join(library_pathname, args.library)
        os.makedirs(library_pathname, exist_ok=True)
        os.chdir(library_pathname)
        url = "blast/db/FASTA/" + args.library + ".gz"
        http_download_file2(
            NCBI_SERVER, [url], save_to=os.path.abspath(os.curdir)
        )
        with gzip.open(args.library + ".gz", mode="rt") as in_file:
            with open(library_filename, "w") as out_file:
                shutil.copyfileobj(in_file, out_file)
        with open(library_filename, "r") as in_file:
            with open("prelim_map.txt", "w") as out_file:
                out_file.write("# prelim_map for " + args.library + "\n")
                scan_fasta_file(
                    in_file,
                    out_file,
                    sequence_to_url="blast/db/FASTA/nt.gz",
                    lenient=True,
                )
    elif args.library in ["UniVec", "UniVec_Core"]:
        if args.protein:
            LOG.error(
                "{:s} is available for nucleotide databases only\n".format(
                    args.library
                )
            )
            sys.exit(1)
        library_pathname = os.path.join(library_pathname, args.library)
        os.makedirs(library_pathname, exist_ok=True)
        os.chdir(library_pathname)
        http_download_file2(
            NCBI_SERVER,
            ["pub/UniVec/" + args.library],
            save_to=os.path.abspath(os.curdir),
        )
        special_taxid = 28384
        LOG.info(
            "Assigning taxonomy ID of {:d} to all sequences\n".format(
                special_taxid
            )
        )
        with open(args.library, "r") as in_file:
            with open("library.fna", "w") as out_file:
                for line in in_file:
                    if line.startswith(">"):
                        line = re.sub(
                            ">",
                            ">kraken:taxid|" + str(special_taxid) + "|",
                            line,
                        )
                    out_file.write(line)
        with open("library.fna", "r") as in_file:
            with open("prelim_map.txt", "w") as out_file:
                out_file.write("# prelim_map for " + args.library + "\n")
                scan_fasta_file(
                    in_file,
                    out_file,
                    sequence_to_url="pub/UniVec/" + args.library,
                )
    else:
        if args.library.upper().startswith("GCF")\
           or args.library.upper().startswith("GCA"):
            download_dataset_by_project(args, "accession",
                                        [args.library.upper()])
        elif args.library.upper().startswith("PRJ"):
            download_dataset_by_project(args, "bioproject",
                                        [args.library.upper()])
        else:
            download_dataset_by_project(args, "taxon", [args.library])

    if not args.no_masking\
       and args.library in ["UniVec", "UniVec_Core", "nt", "nr"]:
        mask_files(
            [library_filename],
            library_filename + ".masked",
            args.masker_threads,
            args.protein
        )
        shutil.move(library_filename + ".masked", library_filename)
    LOG.info("Added {:s} to {:s}\n".format(args.library, args.db))


def get_abs_path(filename):
    return os.path.abspath(filename)


def is_compressed(filename):
    bzip_magic = b"\x42\x5A\x68"
    gzip_magic = b"\x1F\x8B"
    xz_magic = b"\xFD\x37\x7A\x58\x5A\x00"

    nbytes = len(xz_magic)
    with open(filename, "rb") as f:
        data = f.read(nbytes)
        if data.startswith((bzip_magic, gzip_magic, xz_magic)):
            return True
        return False


def get_reader(filename):
    bzip_magic = b"\x42\x5A\x68"
    gzip_magic = b"\x1F\x8B"
    xz_magic = b"\xFD\x37\x7A\x58\x5A\x00"

    nbytes = len(xz_magic)
    with open(filename, "rb") as f:
        data = f.read(nbytes)
        if data.startswith(bzip_magic):
            return bz2.open
        elif data.startswith(gzip_magic):
            return gzip.open
        elif data.startswith(xz_magic):
            return lzma.open
        else:
            return open


def read_from_files(filename1, filename2=None):
    reader1 = get_reader(filename1)
    reader2 = None

    if filename2 is not None:
        reader2 = get_reader(filename2)

    if reader2 is None:
        with reader1(filename1, "rb") as f:
            for seq in f:
                yield seq
    else:
        with reader1(filename1, "rb") as f1, reader2(filename2, "rb") as f2:
            for seq1, seq2 in itertools.zip_longest(f1, f2):
                if seq1 is None:
                    LOG.error(
                        "{} contains more sequences than {}".format(
                            filename1, filename2
                        )
                    )
                    sys.exit(1)
                if seq2 is None:
                    LOG.error(
                        "{} contains more sequences than {}".format(
                            filename2, filename1
                        )
                    )
                    sys.exit(1)
                yield (seq1, seq2)


def write_to_fifo(filenames, fifo1=None, fifo2=None):
    if fifo2 is not None:
        with open(fifo1, "wb") as file1, open(fifo2, "wb") as file2:
            for fn1, fn2 in zip(filenames[0::2], filenames[1::2]):
                for seq1, seq2 in read_from_files(fn1, fn2):
                    file1.write(seq1)
                    file2.write(seq2)
    else:
        with open(fifo1, "wb") as file1:
            for fn in filenames:
                for seq in read_from_files(fn):
                    file1.write(seq)


def check_seqidmap():
    LOG.info(
        "Checking if there are invalid taxid in seqid2taxid.map. "
        "These taxids will be logged if found and removed from the file\n"
    )
    taxonomy_nodes = {}
    with open(os.path.join("taxonomy", "nodes.dmp"), "r") as fin:
        for entry in fin:
            taxid, parent_taxid = entry.split("\t|")[:2]
            taxonomy_nodes[taxid.strip()] = parent_taxid.strip()

    with open("seqid2taxid.map.new", "w") as fout:
        with open("seqid2taxid.map", "r") as fin:
            for line in fin:
                seqid, taxid = line.split("\t")
                taxid = taxid.strip()
                if taxid in taxonomy_nodes:
                    fout.write(line)
                else:
                    LOG.warning(
                        "There is no entry for taxid, '{}', contained in, {},"
                        "in nodes.dmp. Please contact NCBI about this\n"
                        .format(seqid, taxid)
                    )
    shutil.move("seqid2taxid.map.new", "seqid2taxid.map")


def suffix_to_multiplier(suffix):
    name_to_size = {
        'byte': 1,
        'kebibyte': 2 ** 10,
        'mebibyte': 2 ** 20,
        'gebibyte': 2 ** 30,
        'tebibyte': 2 ** 40,
        'kilobyte': 10 ** 3,
        'megabyte': 10 ** 6,
        'gigabyte': 10 ** 9,
        'terabyte': 10 ** 12

    }

    unit_to_size = {
        'B': 1,
        'KiB': 2 ** 10,
        'KB':  10 ** 3,
        'MiB': 2 ** 20,
        'MB':  10 ** 6,
        'GiB': 2 ** 30,
        'GB':  10 ** 9,
        'TiB': 2 ** 40,
        'TB':  10 ** 12,
    }

    original_suffix = suffix

    if suffix in unit_to_size:
        return unit_to_size[suffix]

    if suffix.lower().endswith("s"):
        suffix = suffix.lower()[:-1]

    if suffix in name_to_size:
        return name_to_size[suffix]

    LOG.error("Unable to convert {} to a storage unit\n".format(original_suffix))
    sys.exit(1)



def parse_db_size(input):
    if input.isdigit():
        return int(input)

    input = input.replace(" ", "")
    number = "".join(itertools.takewhile(str.isnumeric, input))
    suffix = "".join(itertools.takewhile(str.isalpha, input[len(number):]))
    number = int(number)
    multiplier = suffix_to_multiplier(suffix)

    return number * multiplier


def build_kraken2_db(args):
    if not os.path.isdir(get_abs_path(args.db)):
        LOG.error('Cannot find Kraken 2 database: "{:s}\n'.format(args.db))
        sys.exit(1)
    os.chdir(args.db)
    if not os.path.isdir("taxonomy"):
        LOG.error("Cannot find taxonomy subdirectory in database\n")
        sys.exit(1)
    if not os.path.isdir("library"):
        LOG.error("Cannot find library subdirectory in database\n")
        sys.exit(1)

    prelim_map_filepaths = []
    prelim_map_mtime = 0
    if os.path.isdir("library"):
        glob_path = os.path.join("library", "*")
        prelim_map_filepaths = glob.glob(
            os.path.join(glob_path, "prelim_map*.txt")
        )
        for prelim_map_filepath in prelim_map_filepaths:
            mtime = os.path.getmtime(prelim_map_filepath)
            if mtime > prelim_map_mtime:
                prelim_map_mtime = mtime

    if os.path.exists("seqid2taxid.map") and \
       os.path.getmtime("seqid2taxid.map") > prelim_map_mtime:
        LOG.info(
            "A seqid2taxid.map already present and newer"
            " than any of the prelim_map.txt files, skipping\n"
        )
    else:
        LOG.info("Concatenating prelim_map.txt files\n")
        with open("prelim_map.txt", "w") as out_file:
            for prelim_map_filepath in prelim_map_filepaths:
                with open(prelim_map_filepath, "r") as in_file:
                    shutil.copyfileobj(in_file, out_file)
        if os.path.getsize("prelim_map.txt") == 0:
            os.remove("prelim_map.txt")
            LOG.error(
                "No preliminary seqid/taxid mapping files found, aborting\n"
            )
            sys.exit(1)
        LOG.info("Finished concatenating prelim_map.txt files\n")
        total_sequences = 0
        LOG.info("Creating sequence ID to taxonomy ID map\n")
        with open("prelim_map.txt", "r") as in_file:
            with open("seqid2taxid.map.tmp", "w") as seqid2taxid_file:
                with open("accmap.tmp", "w") as accmap_file:
                    for line in in_file:
                        if line.startswith("#"):
                            continue
                        line = line.strip()
                        new_line = "\t".join(line.split("\t")[1:3]) + "\n"
                        if line.startswith("TAXID"):
                            seqid2taxid_file.write(new_line)
                        elif line.startswith("ACCNUM"):
                            accmap_file.write(new_line)
                        total_sequences += 1
        if os.path.getsize("accmap.tmp") > 0:
            accession2taxid_filenames = glob.glob("taxonomy/*.accession2taxid")
            if accession2taxid_filenames:
                lookup_accession_numbers(
                    "accmap.tmp",
                    "seqid2taxid.map.tmp",
                    *accession2taxid_filenames
                )
            else:
                LOG.error(
                    "Accession to taxid map files are required to"
                    " build this database.\n"
                )
                LOG.error(
                    "Run k2 download-taxonomy --db {:s} again".format(args.db)
                )
                sys.exit(1)
        os.remove("accmap.tmp")
        move("seqid2taxid.map.tmp", "seqid2taxid.map")
        LOG.info("Created sequence ID to taxonomy ID map\n")
        check_seqidmap()

    estimate_capacity_binary = find_kraken2_binary("estimate_capacity")
    argv = [estimate_capacity_binary, "-S", construct_seed_template(args)]
    if args.protein:
        argv.append("-X")
    wrapper_args_to_binary_args(
        args, argv, get_binary_options(estimate_capacity_binary)
    )
    fasta_filenames = glob.glob(
        os.path.join("library", os.path.join("*", "*.f[an]a")),
        recursive=False
    )
    estimate = ""
    if os.path.exists("estimated_capacity"):
        estimated_capacity_mtime = \
            os.path.getmtime("estimated_capacity")
        seqid_to_taxid_map_mtime = os.path.getmtime("seqid2taxid.map")
        if estimated_capacity_mtime > seqid_to_taxid_map_mtime:
            LOG.info(
                "An estimated_capacity file exists and is newer "
                "than seqid2taxid.map , reading the estimated "
                "capacity from estimated_capacity file.\n"
            )
            with open("estimated_capacity", "r") as in_file:
                estimate = in_file.read()

    if estimate == "":
        if not dwk2():
            argv.extend(fasta_filenames)
        LOG.info("Running: " + " ".join(argv) + "\n")
        proc = subprocess.Popen(
            argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE
        )
        if dwk2():
            for filename in fasta_filenames:
                with open(filename, "rb") as in_file:
                    while True:
                        data = in_file.read(8192)
                        if not data:
                            break
                        proc.stdin.write(data)
        estimate = proc.communicate()[0].decode()
        proc.stdin.close()
        with open("estimated_capacity", "w") as out_file:
            out_file.write(estimate)

    required_capacity = (int(estimate.strip()) + 8192) / args.load_factor
    LOG.info(
        "Estimated hash table requirement: {:s}\n".format(
            format_bytes(required_capacity * 4)
        )
    )
    if args.max_db_size:
        args.max_db_size = parse_db_size(args.max_db_size)
        if args.max_db_size < required_capacity * 4:
            args.max_db_size = int(args.max_db_size / 4)
            LOG.info(
                "Maximum hash table size of {}, specified and is"
                " lower than the calculated estimated capacity of {}\n"
                .format(
                    format_bytes(args.max_db_size * 4),
                    format_bytes(required_capacity * 4))
            )
    if os.path.isfile("hash.k2d"):
        LOG.info("Hash table already present, skipping build\n")
    else:
        LOG.info("Starting database build\n")
        build_db_bin = find_kraken2_binary("build_db")
        argv = [
            build_db_bin,
            "-H",
            "hash.k2d.tmp",
            "-t",
            "taxo.k2d.tmp",
            "-o",
            "opts.k2d.tmp",
            "-n",
            "taxonomy",
            "-m",
            "seqid2taxid.map",
            "-c",
            str(required_capacity),
            "-S",
            construct_seed_template(args),
        ]
        if args.protein:
            argv.append("-X")
        wrapper_args_to_binary_args(
            args, argv, get_binary_options(build_db_bin)
        )

        LOG.info("Running: " + " ".join(argv) + "\n")
        total_sequences = 0
        with open("seqid2taxid.map") as fin:
            for line in fin:
                total_sequences += 1

        m_err, s_err = pty.openpty()
        proc = subprocess.Popen(
            argv, stdin=subprocess.PIPE,
            stdout=s_err,
            stderr=s_err,
        )
        thread = threading.Thread(
            target=read_from_stderr,
            args=(m_err, total_sequences)
        )
        thread.start()
        for filename in fasta_filenames:
            with open(filename, "rb") as in_file:
                while True:
                    data = in_file.read(8192)
                    if not data:
                        break
                    proc.stdin.write(data)

        proc.stdin.close()
        if proc.wait() != 0:
            LOG.error("Encountered error while building database\n")
            sys.exit(1)
        thread.join()
        os.close(m_err)
        os.close(s_err)

        move("hash.k2d.tmp", "hash.k2d")
        move("taxo.k2d.tmp", "taxo.k2d")
        move("opts.k2d.tmp", "opts.k2d")
        LOG.info("Finished building database\n")


def decompress_with_zlib(filename):
    inflator = zlib.decompressobj(15 + 32)
    with open(filename, "rb") as infile:
        while True:
            data = infile.read(8196)
            if not data:
                break
            inflator.decompress(data)


def read_from_stderr(fd, total_sequences):
    pb = ProgressBar(total_sequences)
    buffer = b""
    processing = True
    while processing:
        data = os.read(fd, 1024)
        if len(data) == 0:
            processing = False
        data = buffer + data
        buffer = b""
        for line in data.splitlines(True):
            fields = line.split()
            if line.startswith(b"Processed") and len(fields) > 1:
                buffer = b""
                progress = int(fields[1])
                pb.progress(progress)
                eol = '\r'
                if pb.current == total_sequences:
                    processing = False
                    eol = '\n'
                LOG.debug(
                    "Processed:" +
                    pb.get_bar() + " {}/{}{}"
                    .format(progress, total_sequences, eol)
                )
            elif line.endswith(b"\n"):
                LOG.debug(line.decode())
            else:
                buffer = line


# Parses RDP sequence data to create Kraken taxonomy
# and sequence ID -> taxonomy ID mapping
def build_rdp_taxonomy(f):
    seqid_map = {}
    seen_it = {}
    child_data = {"root;no rank": {}}

    for line in f:
        if not line.startswith(">"):
            continue
        line = line.strip()
        seq_label, taxonomy_string = line.split("\t")
        seqid = seq_label.split(" ")[0]
        taxonomy_string = re.sub(
            "^Lineage=Root;rootrank;", "root;no rank;", taxonomy_string
        )
        taxonomy_string = re.sub(";$", ";no rank", taxonomy_string)
        seqid_map[seqid] = taxonomy_string
        seen_it.setdefault(taxonomy_string, 0)
        seen_it[taxonomy_string] += 1
        if seen_it[taxonomy_string] > 1:
            continue
        while True:
            match = re.search("(;[^;]+;[^;]+)$", taxonomy_string)
            if match is None:
                break
            level = match.group(1)
            taxonomy_string = re.sub(";[^;]+;[^;]+$", "", taxonomy_string)
            key = taxonomy_string + level
            child_data.setdefault(taxonomy_string, {})
            seen_it.setdefault(taxonomy_string, 0)
            child_data[taxonomy_string].setdefault(key, 0)
            child_data[taxonomy_string][key] += 1
            seen_it[taxonomy_string] += 1
            if seen_it[taxonomy_string] > 1:
                break
    id_map = {}
    next_node_id = 1
    with open("names.dmp", "w") as names_file:
        with open("nodes.dmp", "w") as nodes_file:
            bfs_queue = [["root;no rank", 1]]
            while len(bfs_queue) > 0:
                node, parent_id = bfs_queue.pop()
                match = re.search("([^;]+);([^;]+)$", node)
                if match is None:
                    LOG.error(
                        'BFS processing encountered formatting eror, "{:s}"\n'
                        .format(node)
                    )
                    sys.exit(1)
                display_name, rank = match.group(1), match.group(2)
                if rank == "domain":
                    rank = "superkingdom"
                node_id, next_node_id = next_node_id, next_node_id + 1
                id_map[node] = node_id
                names_file.write(
                    "{:d}\t|\t{:s}\t|\t-\t|\tscientific name\t|\n".format(
                        node_id, display_name
                    )
                )
                nodes_file.write(
                    "{:d}\t|\t{:d}\t|\t{:s}\t|\t-\t|\n".format(
                        node_id, parent_id, rank
                    )
                )
                children = (
                    sorted([key for key in child_data[node]])
                    if node in child_data
                    else []
                )
                for node in children:
                    bfs_queue.insert(0, [node, node_id])
    with open("seqid2taxid.map", "w") as f:
        for seqid in sorted([key for key in seqid_map]):
            taxid = id_map[seqid_map[seqid]]
            f.write("{:s}\t{:d}\n".format(seqid, taxid))


# Build the standard Kraken database
def build_standard_database(args):
    download_taxonomy(args)
    for library in [
        "archaea",
        "bacteria",
        "viral",
        "plasmid",
        "human",
        "UniVec_Core",
    ]:
        if library == "UniVec_Core" and args.protein:
            continue
        args.library = library
        download_genomic_library(args)
    build_kraken2_db(args)


# Parses Silva taxonomy file to create Kraken taxonomy
def build_silva_taxonomy(in_file):
    id_map = {"root": 1}
    with open("names.dmp", "w") as names_file:
        with open("nodes.dmp", "w") as nodes_file:
            names_file.write("1\t|\troot\t|\t-\t|\tscientific name\t|\n")
            nodes_file.write("1\t|\t1\t|\tno rank\t|\t-\t|\n")
            for line in in_file:
                line = line.strip()
                taxonomy_string, node_id, rank = line.split("\t")[:3]
                id_map[taxonomy_string] = node_id
                match = re.search("^(.+;|)([^;]+);$", taxonomy_string)
                if match:
                    parent_name = match.group(1)
                    display_name = match.group(2)
                    if parent_name == "":
                        parent_name = "root"
                    parent_id = id_map[parent_name] or None
                    if not parent_id:
                        LOG.error('orphan error: "{:s}"\n'.format(line))
                        sys.exit(1)
                    if rank == "domain":
                        rank = "superkingdom"
                    names_file.write(
                        "{:s}\t|\t{:s}\t|\t-\t|\tscientific name\t|\n".format(
                            node_id, display_name
                        )
                    )
                    nodes_file.write(
                        "{:s}\t|\t{:s}\t|\t{:s}\t|\t-\t|\n".format(
                            node_id, str(parent_id), rank
                        )
                    )
                else:
                    LOG.error('strange input: "{:s}"\n'.format(line))
                    sys.exit(1)


# Build a 16S database from Silva data
def build_16S_silva(args):
    args.db = os.path.abspath(args.db)
    os.makedirs(args.db, exist_ok=True)
    os.chdir(args.db)
    for directory in ["data", "taxonomy", "library"]:
        os.makedirs(directory, exist_ok=True)
    os.chdir("data")
    remote_directory = "/release_138_2/Exports"
    fasta_filename = "SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"
    taxonomy_prefix = "tax_slv_ssu_138.2"
    ftp = FTP(SILVA_SERVER)
    ftp.download(remote_directory, fasta_filename)
    ftp.download(
        remote_directory + "/taxonomy", taxonomy_prefix + ".acc_taxid.gz"
    )
    decompress_files([taxonomy_prefix + ".acc_taxid.gz"])
    ftp.download(remote_directory + "/taxonomy", taxonomy_prefix + ".txt.gz")
    with gzip.open(taxonomy_prefix + ".txt.gz", "rt") as f:
        build_silva_taxonomy(f)
    os.chdir(os.path.pardir)
    move(os.path.join("data", "names.dmp"), "taxonomy")
    move(os.path.join("data", "nodes.dmp"), "taxonomy")
    move(
        os.path.join("data", taxonomy_prefix + ".acc_taxid"), "seqid2taxid.map"
    )
    with gzip.open(os.path.join("data", fasta_filename), "rt") as in_file:
        os.chdir("library")
        os.makedirs("silva", exist_ok=True)
        os.chdir("silva")
        with open("library.fna", "w") as out_file:
            for line in in_file:
                if not line.startswith(">"):
                    line = line.replace("U", "T")
                out_file.write(line)
    if not args.no_masking:
        filename = "library.fna"
        mask_files(
            [filename], filename + ".masked", args.threads
        )
        shutil.move(filename + ".masked", filename)

    os.chdir(args.db)
    build_kraken2_db(args)


# Parses Greengenes taxonomy file to create Kraken taxonomy
# and sequence ID -> taxonomy ID mapping
# Input: gg_13_5_taxonomy.txt
def build_gg_taxonomy(in_file):
    rank_codes = {
        "k": "superkingdom",
        "p": "phylum",
        "c": "class",
        "o": "order",
        "f": "family",
        "g": "genus",
        "s": "species",
    }
    seqid_map = {}
    seen_it = {}
    child_data = {"root": {}}
    for line in in_file:
        line = line.strip()
        seqid, taxonomy_string = line.split("\t")
        taxonomy_string = re.sub("(; [a-z]__)+$", "", taxonomy_string)
        seqid_map[seqid] = taxonomy_string
        seen_it.setdefault(taxonomy_string, 0)
        seen_it[taxonomy_string] += 1
        if seen_it[taxonomy_string] > 1:
            continue
        while True:
            match = re.search("(; [a-z]__[^;]+$)", taxonomy_string)
            if not match:
                break
            level = match.group(1)
            taxonomy_string = re.sub("(; [a-z]__[^;]+$)", "", taxonomy_string)
            child_data.setdefault(taxonomy_string, {})
            key = taxonomy_string + level
            seen_it.setdefault(taxonomy_string, 0)
            child_data[taxonomy_string].setdefault(key, 0)
            child_data[taxonomy_string][key] += 1
            seen_it[taxonomy_string] += 1
            if seen_it[taxonomy_string] > 1:
                break
        if seen_it[taxonomy_string] == 1:
            child_data["root"].setdefault(taxonomy_string, 0)
            child_data["root"][taxonomy_string] += 1
    id_map = {}
    next_node_id = 1
    with open("names.dmp", "w") as names_file:
        with open("nodes.dmp", "w") as nodes_file:
            bfs_queue = [["root", 1]]
            while len(bfs_queue) > 0:
                node, parent_id = bfs_queue.pop()
                display_name = node
                rank = None
                match = re.search("g__([^;]+); s__([^;]+)$", node)
                if match:
                    genus, species = match.group(1), match.group(2)
                    rank = "species"
                    if re.search(" endosymbiont ", species):
                        display_name = species
                    else:
                        display_name = genus + " " + species
                else:
                    match = re.search("([a-z])__([^;]+)$", node)
                    if match:
                        rank = rank_codes[match.group(1)]
                        display_name = match.group(2)
                rank = rank or "no rank"
                node_id, next_node_id = next_node_id, next_node_id + 1
                id_map[node] = node_id
                names_file.write(
                    "{:d}\t|\t{:s}\t|\t-\t|\tscientific name\t|\n".format(
                        node_id, display_name
                    )
                )
                nodes_file.write(
                    "{:d}\t|\t{:d}\t|\t{:s}\t|\t-\t|\n".format(
                        node_id, parent_id, rank
                    )
                )
                children = (
                    sorted([key for key in child_data[node]])
                    if node in child_data
                    else []
                )
                for node in children:
                    bfs_queue.insert(0, [node, node_id])
    with open("seqid2taxid.map", "w") as f:
        for seqid in sorted([key for key in seqid_map], key=int):
            taxid = id_map[seqid_map[seqid]]
            f.write("{:s}\t{:d}\n".format(seqid, taxid))


# Build a 16S database from Greengenes data
def build_16S_gg(args):
    args.db = os.path.abspath(args.db)
    os.makedirs(args.db, exist_ok=True)
    gg_version = "gg_13_5"
    remote_directory = "/greengenes_release/" + gg_version
    os.chdir(args.db)
    for directory in ["data", "taxonomy", "library"]:
        os.makedirs(directory, exist_ok=True)
    os.chdir("data")
    ftp = FTP(GREENGENES_SERVER)
    ftp.download(remote_directory, gg_version + ".fasta.gz")
    decompress_files([gg_version + ".fasta.gz"])
    ftp.download(remote_directory, gg_version + "_taxonomy.txt.gz")
    decompress_files([gg_version + "_taxonomy.txt.gz"])
    with open(gg_version + "_taxonomy.txt", "r") as f:
        build_gg_taxonomy(f)
    os.chdir(os.path.abspath(os.path.pardir))
    move(os.path.join("data", "names.dmp"), "taxonomy")
    move(os.path.join("data", "nodes.dmp"), "taxonomy")
    move(os.path.join("data", "seqid2taxid.map"), os.getcwd())
    move(
        os.path.join("data", gg_version + ".fasta"),
        os.path.join("library", "library.fna"),
    )
    os.chdir("library")
    os.makedirs("greengenes", exist_ok=True)
    move("library.fna", "greengenes")
    os.chdir("greengenes")
    if not args.no_masking:
        filename = "library.fna"
        mask_files([filename], filename + ".masked", args.threads)
        move(filename + ".masked", filename)
    os.chdir(args.db)
    build_kraken2_db(args)


# Build a 16S data from RDP data
def build_16S_rdp(args):
    os.makedirs(args.db, exist_ok=True)
    os.chdir(args.db)
    for directory in ["data", "taxonomy", "library"]:
        os.makedirs(directory, exist_ok=True)
    os.chdir("data")
    http_download_file(
        "http://rdp.cme.msu.edu/download/current_Bacteria_unaligned.fa.gz"
    )
    http_download_file(
        "http://rdp.cme.msu.edu/download/current_Archaea_unaligned.fa.gz"
    )
    decompress_files(glob.glob("*gz"))
    for filename in glob.glob("current_*_unaligned.fa"):
        with open(filename, "r") as f:
            build_rdp_taxonomy(f)
    os.chdir(os.pardir)
    move(os.path.join("data", "names.dmp"), "taxonomy")
    move(os.path.join("data", "nodes.dmp"), "taxonomy")
    move(os.path.join("data", "seqid2taxid.map"), os.getcwd())
    for filename in glob.glob(os.path.join("data", "*.fa")):
        new_filename = os.path.basename(re.sub(r"\.fa$", ".fna", filename))
        shutil.move(filename, os.path.join("library", new_filename))
        if not args.no_masking:
            new_filename = os.path.join("library", new_filename)
            mask_files(
                [new_filename], new_filename + ".masked", args.threads
            )
            shutil.move(new_filename + ".masked", new_filename)

    build_kraken2_db(args)


# Reads multi-FASTA input and examines each sequence header.  In quiet
# mode headers are OK if a taxonomy ID is found (as either the entire
# sequence ID or as part of a "kraken:taxid" token), or if something
# looking like a GI or accession number is found.  In normal mode, the
# taxonomy ID will be looked up (if not explicitly specified in the
# sequence ID) and reported if it can be found.  Output is
# tab-delimited lines, with sequence IDs in first column and taxonomy
# IDs in second.


# Sequence IDs with a kraken:taxid token will use that to assign taxonomy
# ID, e.g.:
# >gi|32499|ref|NC_021949.2|kraken:taxid|562|
#
# Sequence IDs that are completely numeric are assumed to be the taxonomy
# ID for that sequence.
#
# Otherwise, an accession number is searched for; if not found, a GI
# number is searched for.  Failure to find any of the above is a fatal error.
# Without `quiet`, a comma-separated file list specified by -A (for both accession
# numbers and GI numbers) is examined; failure to find a
# taxonomy ID that maps to a provided accession/GI number is non-fatal and
# will emit a warning.
#
# With -q, does not print any output, and will die w/ nonzero exit instead
# of warning when unable to find a taxid, accession #, or GI #.
#
def make_seqid_to_taxid_map(
    in_file, quiet, accession_map_filenames=False, library_map_filename=None
):
    target_lists = {}
    for line in in_file:
        match = re.match(r">(\S+)", line)
        if match is None:
            continue
        seqid = match.group(1)
        output = None
        regexes = [
            r"(?:^|\|)kraken:taxid\|(\d+)",
            r"^\d+$",
            r"(?:^|\|)([A-Z]+_?[A-Z0-9]+)(?:\||\b|\.)",
            r"(?:^|\|)gi\|(\d+)",
        ]
        match = None
        index = None
        for i, regex in enumerate(regexes):
            match = re.match(regex, seqid)
            if match:
                index = i
                break
        if index == 0:
            output = seqid + "\t" + match.group(1) + "\n"
        elif index == 1:
            output = seqid + "\t" + seqid + "\n"
        elif index in [2, 3]:
            if not quiet:
                capture = match.group(1)
                target_lists.setdefault(capture, [])
                target_lists[capture].insert(0, seqid)
        else:
            LOG.error(
                "Unable to determine taxonomy ID for sequence {:s}\n".format(
                    seqid
                )
            )
            sys.exit(1)
        if output and not quiet:
            print(output)
    if quiet:
        if len(target_lists) == 0:
            LOG.error("External map required\n")
        sys.exit(0)
    if len(target_lists) == 0:
        sys.exit(0)
    if not accession_map_filenames and library_map_filename is None:
        LOG.error(
            "Found sequence ID without explicit taxonomy ID, but no map used\n"
        )
        sys.exit(1)
    # Remove targets where we've already handled the mapping
    if library_map_filename:
        with open(library_map_filename, "r") as f:
            for line in f:
                line = line.strip()
                seqid, taxid = line.split("\t")
                if seqid in target_lists:
                    print("{:s}\t{:s}\n".format(seqid, taxid))
                    del target_lists[seqid]
    if len(target_lists) == 0:
        sys.exit(0)
    for filename in accession_map_filenames:
        with open(filename, "r") as f:
            f.readline()
            for line in f:
                line = line.strip()
                accession, with_version, taxid, gi = line.split("\t")
                if accession in target_lists:
                    target_list = target_lists[accession]
                    del target_lists[accession]
                    for seqid in target_list:
                        print("{:s}\t{:s}".format(seqid, taxid))
                if gi != "na" and gi in target_lists:
                    target_list = target_lists[gi]
                    del target_lists[gi]
                    for seqid in target_list:
                        print("{:s}\t{:s}\n".format(seqid, taxid))


def classify(args):
    classify_bin = find_kraken2_binary("classify")
    database_path = find_database(args.db)
    if database_path is None:
        LOG.error("{:s} is not a valid database... exiting\n".format(args.db))
        sys.exit(1)
    if "paired" in args and len(args.filenames) % 2 != 0:
        LOG.error("--paired requires an even number of file names\n")
        sys.exit(1)
    if args.confidence < 0 or args.confidence > 1:
        LOG.error(
            "--confidence, {:f}, must be between 0 and 1 inclusive\n".format(
                args.confidence
            )
        )
        sys.exit(1)
    argv = [
        classify_bin,
        "-H",
        os.path.join(database_path, "hash.k2d"),
        "-t",
        os.path.join(database_path, "taxo.k2d"),
        "-o",
        os.path.join(database_path, "opts.k2d"),
    ]
    wrapper_args_to_binary_args(args, argv, get_binary_options(classify_bin))
    if any([is_compressed(filename) for filename in args.filenames]):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            fifo1_pathname = os.path.join(temp_dir_name, "fifo1")
            fifo2_pathname = None
            try:
                os.mkfifo(fifo1_pathname, 0o600)
            except OSError:
                LOG.error(
                    "Unable to create FIFO for processing compressed files\n"
                )
                sys.exit(1)
            if "-P" in argv:
                fifo2_pathname = os.path.join(temp_dir_name, "fifo2")
                try:
                    os.mkfifo(fifo2_pathname, 0o600)
                except OSError:
                    LOG.error(
                        "Unable to create FIFO for processing compressed files\n"
                    )
                    sys.exit(1)
                argv.extend([fifo1_pathname, fifo2_pathname])
            else:
                argv.append(fifo1_pathname)
            thread = threading.Thread(target=subprocess.call, args=(argv,))
            thread.start()
            write_to_fifo(args.filenames, fifo1_pathname, fifo2_pathname)
            thread.join()
    else:
        argv.extend(args.filenames)
        sys.exit(subprocess.call(argv))


def inspect_db(args):
    database_pathname = find_database(args.db)
    if not database_pathname:
        LOG.error("{:s} database does not exist\n".format(args.db))
        sys.exit(1)
    for database_file in ["taxo.k2d", "hash.k2d", "opts.k2d"]:
        if not os.path.isfile(os.path.join(database_pathname, database_file)):
            LOG.error("{:s} does not exist\n".format(database_file))
    dump_table_bin = find_kraken2_binary("dump_table")
    argv = [
        dump_table_bin,
        "-H",
        os.path.join(database_pathname, "hash.k2d"),
        "-t",
        os.path.join(database_pathname, "taxo.k2d"),
        "-o",
        os.path.join(database_pathname, "opts.k2d"),
    ]
    # dump_table does not save the table header to file.
    # This is a workaround helps enables us to capture
    # the entire output.
    output_filename, args.output = args.output, None
    wrapper_args_to_binary_args(args, argv, get_binary_options(dump_table_bin))
    process = subprocess.Popen(
        argv, stdout=subprocess.PIPE
    )

    if output_filename == "-":
        shutil.copyfileobj(process.stdout, sys.stdout)
    else:
        with open(output_filename, "wb") as fout:
            shutil.copyfileobj(process.stdout, fout)
    sys.exit(process.wait())


def format_bytes(size):
    current_suffix = "B"
    for suffix in ["kB", "MB", "GB", "TB", "PB", "EB"]:
        if size >= 1024:
            current_suffix = suffix
            size /= 1024
        else:
            break
    return "{:.2f}{:s}".format(size, current_suffix)


def clean_up(filenames):
    LOG.info("Removing the following files: {}\n".format(filenames))
    # walk the directory tree to get the size of the individual files
    # sum them up to get the usage stat
    space_freed = format_bytes(remove_files(filenames))
    LOG.info(
        "Cleaned up {} of space\n".format(space_freed)
    )


def clean_db(args):
    os.chdir(args.db)
    if args.pattern:
        clean_up(glob.glob(args.pattern, recursive=False))
    else:
        clean_up(
            [
                "data",
                "library",
                "taxonomy",
                "seqid2taxid.map",
                "prelim_map.txt",
            ]
        )


def make_build_parser(subparsers):
    parser = subparsers.add_parser(
        "build",
        help="Build a database from library\
              (requires taxonomy which can be downloading\
              via download-taxonomy subcommand, and at least one library\
              which can be added via the download-library or\
              add-to-library subcommands).",
    )
    parser.add_argument(
        "--db",
        type=str,
        metavar="PATHNAME",
        required=True,
        help="Pathname to database folder where building will take place.",
    )
    group = parser.add_argument_group("special")
    mutex_group = group.add_mutually_exclusive_group()
    mutex_group.add_argument(
        "--standard",
        action="store_true",
        help="Make standard database which includes: archaea,\
               bacteria, human, plasmid, UniVec_Core, and viral."
    )
    mutex_group.add_argument(
        "--special",
        type=str,
        choices=["greengenes", "rdp", "silva", "gtdb"],
        help="Build special database. RDP is currently unavailable\
              as URLs no longer work.",
    )
    group.add_argument(
        "--gtdb-files",
        type=str,
        nargs="+",
        help="A list of files or regex matching the files needed to build\
              the special database."
    )
    group.add_argument(
        "--no-masking",
        action="store_true",
        help="Avoid masking low-complexity sequences prior to\
              building database.",
    )
    group.add_argument(
        "--masker-threads",
        type=int,
        default=4,
        metavar="K2MASK_THREADS",
        help="Number of threads used by k2mask during masking\
              process (default: 4)"
    )
    parser.add_argument(
        "--kmer-len", type=int, metavar="INT", help="K-mer length in bp/aa"
    )
    parser.add_argument(
        "--minimizer-len",
        type=int,
        metavar="INT",
        help="Minimizer length in bp/aa",
    )
    parser.add_argument(
        "--minimizer-spaces",
        type=int,
        metavar="INT",
        help="Number of characters in minimizer that are\
              ignored in comparisons",
    )
    parser.add_argument(
        "--threads",
        type=int,
        metavar="INT",
        default=os.environ.get("KRAKEN2_NUM_THREADS") or 1,
        help="Number of threads",
    )
    parser.add_argument(
        "--load-factor",
        type=float,
        metavar="FLOAT (0,1]",
        default=0.7,
        help="Proportion of the hash table to be populated (default: 0.7)",
    )
    parser.add_argument(
        "--fast-build",
        action="store_true",
        help="Do not require database to be deterministically\
              built when using multiple threads. This is faster, but\
              does introduce variability in minimizer/LCA pairs.",
    )
    parser.add_argument(
        "--max-db-size",
        # type=int,
        metavar="SIZE",
        help="Maximum number of bytes for Kraken 2 hash table;\
              if the estimator determines more would normally be\
              needed, the reference library will be downsampled to fit",
    )
    parser.add_argument(
        "--skip-maps",
        action="store_true",
        help="Avoids downloading accession number to taxid maps",
    )
    parser.add_argument(
        "--protein",
        action="store_true",
        help="Build a protein database for translated search",
    )
    parser.add_argument(
        "--block-size",
        type=int,
        metavar="INT",
        default=16384,
        help="Read block size (default: 16384)",
    )
    parser.add_argument(
        "--sub-block-size",
        type=int,
        metavar="INT",
        default=0,
        help="Read subblock size",
    )
    parser.add_argument(
        "--minimum-bits-for-taxid",
        type=int,
        metavar="INT",
        default=0,
        help="Bit storage requested for taxid",
    )
    parser.add_argument(
        "--log",
        type=str,
        metavar="FILENAME",
        default=None,
        help="Specify a log file (default: stderr)",
    )


def make_download_taxonomy_parser(subparsers):
    parser = subparsers.add_parser(
        "download-taxonomy", help="Download NCBI taxonomic information"
    )
    parser.add_argument(
        "--db",
        type=str,
        metavar="PATHNAME",
        required=True,
        help="Pathname to Kraken2 database",
    )
    # parser.add_argument(
    #     "--source",
    #     type=str,
    #     choices=[
    #         "GTDB",
    #         "NCBI"
    #     ],
    #     default="NCBI",
    #     help="From which database should the files be downloaded"
    # )
    parser.add_argument(
        "--protein",
        action="store_true",
        help="Files being added are for a protein database",
    )
    parser.add_argument(
        "--skip-maps",
        action="store_true",
        help="Avoids downloading accession number to taxid maps",
    )
    parser.add_argument(
        "--log",
        type=str,
        metavar="FILENAME",
        default=None,
        help="Specify a log filename (default: stderr)",
    )


def make_download_library_parser(subparsers):
    parser = subparsers.add_parser(
        "download-library", aliases=["download"],
        help="Download and build a special database"
    )
    parser.add_argument(
        "--db",
        type=str,
        metavar="PATHNAME",
        required=True,
        help="Pathname to Kraken2 database",
    )
    parser.add_argument(
        "--library",
        "--taxid",
        "--project",
        "--accession",
        type=str,
        dest="library",
        required=True,
        # choices=[
        #     "archaea",
        #     "bacteria",
        #     "plasmid",
        #     "plastid",
        #     "viral",
        #     "human",
        #     "invertebrate",
        #     "fungi",
        #     "plant",
        #     "protozoa",
        #     "vertebrate_other"
        #     "vertebrate_mammalian",
        #     "mitochondrion",
        #     "nr",
        #     "nt",
        #     "UniVec",
        #     "UniVec_Core",
        # ],
        help="Name of library to download",
    )
    parser.add_argument(
        "--assembly-source",
        type=str,
        required=False,
        choices=["refseq", "genbank", "all"],
        default="refseq",
        help="Download RefSeq (GCF_) or GenBank (GCA_) genome assemblies\
              or both (default RefSeq)",
    )
    parser.add_argument(
        "--assembly-levels",
        type=str,
        nargs="+",
        choices=["chromosome", "complete_genome", "scaffold", "contig"],
        default=["chromosome", "complete_genome"],
        help="Only return genome assemblies that have one of the specified\
              assembly levels (default chromosome and complete genome)"
    )
    parser.add_argument(
        "--has-annotation",
        action="store_true",
        help="Return only annotated genome assemblies (default false)"
    )
    parser.add_argument(
        "--protein",
        action="store_true",
        help="Files being added are for a protein database",
    )
    parser.add_argument(
        "--log",
        type=str,
        metavar="FILENAME",
        default=None,
        help="Specify a log filename (default stderr)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        metavar="THREADS",
        default=1,
        help="The number of threads/processes k2 uses when downloading\
              and processing library files.",
    )
    masking_parser = parser.add_mutually_exclusive_group()
    masking_parser.add_argument(
        "--no-masking",
        action="store_true",
        help="Avoid asking low-complexity sequences prior to\
              building; masking requires k2mask or segmasker to be\
              installed",
    )
    masking_parser.add_argument(
        "--masker-threads",
        type=int,
        default=4,
        metavar="K2MASK_THREADS",
        help="Number of threads used by k2mask during masking\
              process (default: 4)"
    )


def make_add_to_library_parser(subparsers):
    parser = subparsers.add_parser(
        "add-to-library", help="Add file(s) to library"
    )
    parser.add_argument(
        "--db",
        type=str,
        metavar="PATHNAME",
        required=True,
        help="Pathname to Kraken2 database",
    )
    parser.add_argument(
        "--threads",
        type=int,
        metavar="THREADS",
        default=1,
        help="The number of threads/processes k2 uses when\
              adding library files."
    )
    parser.add_argument(
        "--file",
        "--files",
        type=str,
        nargs="+",
        required=True,
        dest="files",
        help="""Pathname or patterns of file(s) to be added to library.
                Supported pattern are as follows:
              ? - A question-mark is a pattern that shall match any
                  character.
              * - An asterisk is a pattern that shall match multiple
                  characters.
              [ - The open bracket shall introduce a pattern bracket
                  expression.
             ** - will match any files and zero or more directories,
                  subdirectories and symbolic links to directories.
        """,
    )
    parser.add_argument(
        "--protein",
        action="store_true",
        help="Files being added are for a protein database",
    )
    parser.add_argument(
        "--log",
        type=str,
        metavar="FILENAME",
        default=None,
        help="Specify a log filename (default: stderr)",
    )
    # parser.add_argument(
    #     "--skip-md5",
    #     action="store_true",
    #     help="K2 will by default perform an MD5 check to determine whether\
    #          a file has already been added. This option will allow the user\
    #          to skip this and instead simply compare filenames."
    # )
    masking_parser = parser.add_mutually_exclusive_group()
    masking_parser.add_argument(
        "--no-masking",
        action="store_true",
        help="Avoid asking low-complexity sequences prior to\
              building; masking requires k2mask or segmasker to be\
              installed",
    )
    masking_parser.add_argument(
        "--masker-threads",
        type=int,
        metavar="K2MASK_THREADS",
        default=4,
        help="Number of threads used by k2mask during masking process\
              (default: 4)"
    )


def make_classify_parser(subparsers):
    parser = subparsers.add_parser(
        "classify", help="Classify a set of sequences"
    )
    parser.add_argument(
        "--db",
        type=str,
        metavar="PATHNAME",
        required=True,
        help="Pathname to Kraken2 database.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        metavar="INT",
        default=os.environ.get("KRAKEN2_NUM_THREADS") or 1,
        help="Number of threads",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        default=argparse.SUPPRESS,
        help="Quick operation (use first hit or hits)",
    )
    parser.add_argument(
        "--unclassified-out",
        type=str,
        default=argparse.SUPPRESS,
        metavar="FILENAME",
        help="Print unclassified sequences to filename",
    )
    parser.add_argument(
        "--classified-out",
        type=str,
        metavar="FILENAME",
        default=argparse.SUPPRESS,
        help="Print classified sequences to filename",
    )
    parser.add_argument(
        "--output",
        type=str,
        metavar="FILENAME",
        default=argparse.SUPPRESS,
        help='Print output to file (default: stdout) "-" will \
              suppress normal output',
    )
    parser.add_argument(
        "--confidence",
        type=float,
        default=0.0,
        help="confidence score threshold (default: 0.0); must be in [0,1]",
    )
    parser.add_argument(
        "--minimum-base-quality",
        type=int,
        metavar="INT",
        default=0,
        help="Minimum base quality used in classification",
    )
    parser.add_argument(
        "--report",
        type=str,
        default=argparse.SUPPRESS,
        help="Print a report with aggregate counts/clade to file",
    )
    parser.add_argument(
        "--use-mpa-style",
        action="store_true",
        default=argparse.SUPPRESS,
        help="With --report, format report output like Kraken 1's\
              kraken-mpa-report",
    )
    parser.add_argument(
        "--report-zero-counts",
        action="store_true",
        default=argparse.SUPPRESS,
        help="With --report, report counts for ALL taxa, even if\
              counts are zero",
    )
    parser.add_argument(
        "--report-minimizer-data",
        action="store_true",
        default=argparse.SUPPRESS,
        help="With --report, report minimizer and distinct minimizer\
              count information in addition to normal Kraken report",
    )
    parser.add_argument(
        "--memory-mapping",
        action="store_true",
        default=argparse.SUPPRESS,
        help="Avoids loading entire database into RAM",
    )
    paired_group = parser.add_mutually_exclusive_group()
    paired_group.add_argument(
        "--paired",
        action="store_true",
        default=argparse.SUPPRESS,
        help="The filenames provided have paired-end reads",
    )
    paired_group.add_argument(
        "--interleaved",
        action="store_true",
        default=argparse.SUPPRESS,
        help="The filenames provided have paired-end reads",
    )
    parser.add_argument(
        "--use-names",
        action="store_true",
        default=argparse.SUPPRESS,
        help="Print scientific names instead of just taxids",
    )
    parser.add_argument(
        "--minimum-hit-groups",
        type=int,
        metavar="INT",
        default=2,
        help="Minimum number of hit groups (overlapping k-mers\
              sharing the same minimizer) needed to make a call\
              (default 2)",
    )
    parser.add_argument(
        "--log",
        type=str,
        metavar="FILENAME",
        default=None,
        help="Specify a log filename (default: stderr)",
    )
    parser.add_argument(
        "filenames",
        nargs="+",
        type=str,
        help="Filenames to be classified, supports bz2, gzip, and xz",
    )


def make_inspect_parser(subparsers):
    parser = subparsers.add_parser("inspect", help="Inspect Kraken 2 database")
    parser.add_argument(
        "--db",
        type=str,
        metavar="PATHNAME",
        required=True,
        help="Pathname to Kraken2 database",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=os.environ.get("KRAKEN2_NUM_THREADS") or 1,
        help="Number of threads",
    )
    parser.add_argument(
        "--skip-counts",
        action="store_true",
        help="Only print database summary statistics",
    )
    parser.add_argument(
        "--use-mpa-style",
        action="store_true",
        help="Format output like Kraken 1's kraken-mpa-report",
    )
    parser.add_argument(
        "--report-zero-counts",
        action="store_true",
        help="Report counts for ALL taxa, even if counts are zero",
    )
    parser.add_argument(
        "--log",
        type=str,
        metavar="FILENAME",
        default=None,
        help="Specify a log filename (default: stderr)",
    )
    parser.add_argument(
        "--output", "--out",
        type=str,
        metavar="FILENAME",
        default="-",
        help="Write inspect output to FILENAME (default: stdout)"
    )
    parser.add_argument(
        "--memory-mapping",
        action="store_true",
        default=argparse.SUPPRESS,
        help="Avoids loading entire database into RAM",
    )


def make_clean_parser(subparsers):
    parser = subparsers.add_parser(
        "clean", help="Removes unwanted files from database"
    )
    parser.add_argument(
        "--db",
        type=str,
        metavar="PATHNAME",
        required=True,
        help="Pathname to Kraken2 database",
    )
    parser.add_argument(
        "--log",
        type=str,
        metavar="FILENAME",
        default=None,
        help="Specify a log filename (default: stderr)",
    )
    parser.add_argument(
        "--pattern",
        type=str,
        metavar="SHELL_REGEX",
        default=None,
        help="""Files that match this regular expression will be deleted.
              ? - A question-mark is a pattern that shall match any
                  character.
              * - An asterisk is a pattern that shall match multiple
                  characters.
              [ - The open bracket shall introduce a pattern bracket
                  expression.
             ** - will match any files and zero or more directories,
                  subdirectories and symbolic links to directories.
        """
    )


class HelpAction(argparse._HelpAction):
    def __call__(self, parser, namespace, values, option_string=None):
        parser.print_help()
        subparsers = parser._actions[1].choices
        for action, arg_parser in subparsers.items():
            sys.stderr.write("\n\n" + action + "\n" + "-" * len(action) + "\n")
            arg_parser.print_help()
        sys.exit(0)


def make_cmdline_parser():
    parser = argparse.ArgumentParser("k2", add_help=False)
    parser.add_argument("-h", "--help", action=HelpAction)
    subparsers = parser.add_subparsers()
    make_add_to_library_parser(subparsers)
    make_download_library_parser(subparsers)
    make_download_taxonomy_parser(subparsers)
    make_build_parser(subparsers)
    make_classify_parser(subparsers)
    make_inspect_parser(subparsers)
    make_clean_parser(subparsers)
    return parser


def setup_logger(filename=None):
    logging.StreamHandler.terminator = ""
    logger = logging.getLogger("kraken2")
    if filename:
        logger.setLevel(logging.INFO)
        handler = logging.FileHandler(filename)
        formatter = logging.Formatter("[%(levelname)s - %(asctime)s]: %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    else:
        logger.setLevel(logging.DEBUG)
        handler = logging.StreamHandler()
        formatter = logging.Formatter("[%(levelname)s - %(asctime)s]: %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger


def k2_main():
    global SCRIPT_PATHNAME
    global LOG

    SCRIPT_PATHNAME = os.path.realpath(inspect.getsourcefile(k2_main))

    parser = make_cmdline_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(sys.argv[1:])
    LOG = setup_logger(args.log)
    task = sys.argv[1]
    if task not in ["classify", "inspect"]:
        args.db = os.path.abspath(args.db)
    if task == "download-taxonomy":
        download_taxonomy(args)
    elif task == "classify":
        classify(args)
    elif task == "download-library":
        download_genomic_library(args)
    elif task == "add-to-library":
        add_to_library(args)
    elif task == "inspect":
        inspect_db(args)
    elif task == "clean":
        clean_db(args)
    elif task == "build":
        # Protein defaults
        default_aa_minimizer_length = 12
        default_aa_kmer_length = 15
        default_aa_minimizer_spaces = 0
        # Nucleotide defaults
        default_nt_minimizer_length = 31
        default_nt_kmer_length = 35
        default_nt_minimizer_spaces = 7

        if args.sub_block_size == 0:
            args.sub_block_size = math.ceil(args.block_size / args.threads)
        if not args.kmer_len:
            args.kmer_len = (
                default_aa_kmer_length
                if args.protein
                else default_nt_kmer_length
            )
        if not args.minimizer_len:
            args.minimizer_len = (
                default_aa_minimizer_length
                if args.protein
                else default_nt_minimizer_length
            )
        if not args.minimizer_spaces:
            args.minimizer_spaces = (
                default_aa_minimizer_spaces
                if args.protein
                else default_nt_minimizer_spaces
            )
        if args.minimizer_len > args.kmer_len:
            LOG.error(
                "Minimizer length ({}) must not be greater than kmer "
                "length {}\n".format(args.minimizer_len, args.kmer_len)
            )
            sys.exit(1)
        if args.load_factor <= 0 or args.load_factor > 1:
            LOG.error(
                "Load factor must be greater than 0 but no more than 1\n"
            )
            sys.exit(1)
        if args.minimizer_len <= 0 or args.minimizer_len > 31:
            LOG.error(
                "Minimizer length must be a positive integer "
                "and cannot exceed 31.\n"
            )
            sys.exit(1)
        if args.standard:
            build_standard_database(args)
        elif args.special:
            if args.special == "greengenes":
                build_16S_gg(args)
            elif args.special == "silva":
                build_16S_silva(args)
            elif args.special == "gtdb":
                if not args.gtdb_files:
                    LOG.error("Please specify a list of files or pattern of\
                    the files needed to build a GTDB database.\n")
                    sys.exit(1)
                build_gtdb_database(args)
            else:
                # build_16S_rdp(args)
                LOG.error("RDP database no longer supported.\n")
                sys.exit(1)
        else:
            if args.no_masking:
                LOG.warning(
                    "--no-masking only affects the `--standard` and"
                    "`--special` flags. Its effect will be ignored.\n"
                )
            build_kraken2_db(args)


def sigint(signum, frame):
    sys.exit(1)


if __name__ == "__main__":
    signal.signal(signal.SIGINT, sigint)
    k2_main()
