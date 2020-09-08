Introduction
============

Kraken is a taxonomic sequence classifier that assigns taxonomic
labels to DNA sequences.  Kraken examines the $k$-mers within
a query sequence and uses the information within those $k$-mers
to query a database.  That database maps $k$-mers to the lowest
common ancestor (LCA) of all genomes known to contain a given $k$-mer.

The first version of [Kraken] used a large indexed and sorted list of
$k$-mer/LCA pairs as its database.  While fast, the large memory
requirements posed some problems for users, and so Kraken 2 was
created to provide a solution to those problems.

[Kraken 2] differs from Kraken 1 in several important ways:

1. Only minimizers of the $k$-mers in the query sequences are used
as database queries.  Similarly, only minimizers of the $k$-mers in
the reference sequences in the database's genomic library are stored
in the database.  We will also refer to the minimizers as $\ell$-mers,
where $\ell \leq k$.  All $k$-mers are considered to have the same LCA
as their minimizer's database LCA value.
2. Kraken 2 uses a compact hash table that is a probabilistic data
structure.  This means that occasionally, database queries will fail
by either returning the wrong LCA, or by not resulting in a search
failure when a queried minimizer was never actually stored in the
database.  By incurring the risk of these false positives in the data
structure, Kraken 2 is able to achieve faster speeds and lower memory
requirements.  Users should be aware that database false positive
errors occur in less than 1% of queries, and can be compensated for
by use of confidence scoring thresholds.
3. Kraken 2 has the ability to build a database from amino acid
sequences and perform a translated search of the query sequences
against that database.
4. Kraken 2 utilizes spaced seeds in the storage and querying of
minimizers to improve classification accuracy.
5. Kraken 2 provides support for "special" databases that are
not based on NCBI's taxonomy.  These are currently limited to
three popular 16S databases.

Because Kraken 2 only stores minimizers in its hash table, and $k$ can be
much larger than $\ell$, only a small percentage
of the possible $\ell$-mers in a genomic library are actually deposited in
the database.  This creates a situation similar to the Kraken 1 "MiniKraken"
databases; however, preliminary testing has shown the accuracy of a reduced
Kraken 2 database to be quite similar to the full-sized Kraken 2 database,
while Kraken 1's MiniKraken databases often resulted in a substantial loss
of per-read sensitivity.

If you use Kraken 2 in your own work, please cite either the
[Kraken 2 paper] and/or the original [Kraken paper] as appropriate.  Thank you!

[Kraken]: https://ccb.jhu.edu/software/kraken/
[Kraken 2]: https://ccb.jhu.edu/software/kraken2/
[Kraken paper]: https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-3-r46
[Kraken 2 paper]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0


System Requirements
===================

* **Disk space**: Construction of a Kraken 2 standard database requires
    approximately 100 GB of disk space.  A test on 01 Jan 2018 of the
    default installation showed 42 GB of disk space was used to store
    the genomic library files, 26 GB was used to store the taxonomy
    information from NCBI, and 29 GB was used to store the Kraken 2
    compact hash table.

    Like in Kraken 1, we strongly suggest against using NFS storage
    to store the Kraken 2 database if at all possible.

* **Memory**: To run efficiently, Kraken 2 requires enough free memory
    to hold the database (primarily the hash table) in RAM.  While this
    can be accomplished with a ramdisk, Kraken 2 will by default load
    the database into process-local RAM; the `--memory-mapping` switch
    to `kraken2`  will avoid doing so.  The default database size is 29 GB
    (as of Jan. 2018), and you will need slightly more than that in
    RAM if you want to build the default database.

* **Dependencies**: Kraken 2 currently makes extensive use of Linux
    utilities such as sed, find, and wget.  Many scripts are written
    using the Bash shell, and the main scripts are written using Perl.
    Core programs needed to build the database and run the classifier
    are written in C++11, and need to be compiled using a somewhat
    recent version of g++ that will support C++11.  Multithreading is
    handled using OpenMP.  Downloads of NCBI data are performed by wget
    and rsync.  Most Linux systems will have all of the above listed
    programs and development libraries available either by default or
    via package download.

    Unlike Kraken 1, Kraken 2 does not use an external $k$-mer counter.
    However, by default, Kraken 2 will attempt to use the `dustmasker` or
    `segmasker` programs provided as part of NCBI's BLAST suite to mask
    low-complexity regions (see [Masking of Low-complexity Sequences]).

    **MacOS NOTE:** MacOS and other non-Linux operating systems are *not*
    explicitly supported by the developers, and MacOS users should refer to
    the Kraken-users group for support in installing the appropriate utilities
    to allow for full operation of Kraken 2.  We will attempt to use
    MacOS-compliant code when possible, but development and testing time
    is at a premium and we cannot guarantee that Kraken 2 will install
    and work to its full potential on a default installation of MacOS.

    In particular, we note that the default MacOS X installation of GCC
    does not have support for OpenMP.  Without OpenMP, Kraken 2 is
    limited to single-threaded operation, resulting in slower build and
    classification runtimes.

* **Network connectivity**: Kraken 2's standard database build and download
    commands expect unfettered FTP and rsync access to the NCBI FTP
    server.  If you're working behind a proxy, you may need to set
    certain environment variables (such as `ftp_proxy` or `RSYNC_PROXY`)
    in order to get these commands to work properly.

    Kraken 2's scripts default to using rsync for most downloads; however, you
    may find that your network situation prevents use of rsync. In such cases,
    you can try the `--use-ftp` option to `kraken2-build` to force the
    downloads to occur via FTP.

* **MiniKraken**: At present, users with low-memory computing environments
    can replicate the "MiniKraken" functionality of Kraken 1 in two ways:
    first, by increasing
    the value of $k$ with respect to $\ell$ (using the `--kmer-len` and
    `--minimizer-len` options to `kraken2-build`); and secondly, through
    downsampling of minimizers (from both the database and query sequences)
    using a hash function.  This second option is performed if
    the `--max-db-size` option to `kraken2-build` is used; however, the two
    options are not mutually exclusive.
    In a difference from Kraken 1, Kraken 2 does not require building a full
    database and then shrinking it to obtain a reduced database.


Installation
============

To begin using Kraken 2, you will first need to install it, and then
either download or create a database.

Kraken 2 consists of two main scripts (`kraken2` and `kraken2-build`),
along with several programs and smaller scripts.  As part of the installation
process, all scripts and programs are installed in the same directory.
After installation, you can move the main scripts elsewhere, but moving
the other scripts and programs requires editing the scripts and changing
the `$KRAKEN2_DIR` variables in the main scripts.

Once an install directory is selected, you need to run the following
command in the directory where you extracted the Kraken 2 source:

    ./install_kraken2.sh $KRAKEN2_DIR

(Replace `$KRAKEN2_DIR` above with the directory where you want to install
Kraken 2's programs/scripts.)

The `install_kraken2.sh` script should compile all of Kraken 2's code
and setup your Kraken 2 program directory.  Installation is successful if
you see the message "`Kraken 2 installation complete.`"

Once installation is complete, you may want to copy the main Kraken 2
scripts into a directory found in your `PATH` variable (e.g., "`$HOME/bin`"):

    cp $KRAKEN2_DIR/kraken2{,-build,-inspect} $HOME/bin

After installation, you're ready to either create or download a database.


Kraken 2 Databases
==================

A Kraken 2 database is a directory containing at least 3 files:

* `hash.k2d`: Contains the minimizer to taxon mappings
* `opts.k2d`: Contains information about the options used to build the database
* `taxo.k2d`: Contains taxonomy information used to build the database

None of these three files are in a human-readable format.  Other files
may also be present as part of the database build process, and can, if
desired, be removed after a successful build of the database.

In interacting with Kraken 2, you should not have to directly reference
any of these files, but rather simply provide the name of the directory
in which they are stored.  Kraken 2 allows both the use of a standard
database as well as custom databases; these are described in the
sections [Standard Kraken 2 Database] and [Custom Databases] below,
respectively.


Standard Kraken 2 Database
==========================

To create the standard Kraken 2 database, you can use the following command:

    kraken2-build --standard --db $DBNAME

(Replace "`$DBNAME`" above with your preferred database name/location.
Please note that the database will use approximately 100 GB of
disk space during creation, with the majority of that being reference
sequences or taxonomy mapping information that can be removed after the
build.)

This will download NCBI taxonomic information, as well as the
complete genomes in RefSeq for the bacterial, archaeal, and
viral domains, along with the human genome and a collection of
known vectors (UniVec_Core).  After downloading all this data, the build
process begins; this can be the most time-consuming step.  If you
have multiple processing cores, you can run this process with
multiple threads, e.g.:

    kraken2-build --standard --threads 24 --db $DBNAME

Using 32 threads on an AWS EC2 r4.8xlarge instance with 16 dual-core
hyperthreaded 2.30 GHz CPUs and 244 GB of RAM, the build process took
approximately 35 minutes in Jan. 2018.

The build process itself has two main steps, each of which requires passing
over the contents of the reference library:

  1. **Estimation** of the capacity needed in the Kraken 2 compact hash table.
     This uses a low-memory method to reliably estimate the number of
     minimizers present in the reference library given the selected parameters
     $k$ and $\ell$.
  2. **Population** of the hash table (and conversion of the taxonomy to an
     internal format).  This step is a second pass over the reference library
     to find minimizers and then place them in the database.

(There is one other preliminary step where sequence IDs are mapped to
taxonomy IDs, but this is usually a rather quick process and is mostly handled
during library downloading.)

Unlike Kraken 1's build process, Kraken 2 does not perform checkpointing
after the estimation step.  This is because the estimation step is dependent
on the selected $k$ and $\ell$ values, and if the population step fails, it is
likely because $k$ needs to be increased (reducing the overall memory
requirements).


Classification
==============

To classify a set of sequences, use the `kraken2` command:

    kraken2 --db $DBNAME seqs.fa

Output will be sent to standard output by default.  The files
containing the sequences to be classified should be specified
on the command line.  Sequences can also be provided through
standard input using the special filename `/dev/fd/0`.


The `kraken2` program allows several different options:

* **Multithreading**: Use the `--threads NUM` switch to use multiple
    threads.

* **Quick operation**: Rather than searching all $\ell$-mers in a sequence,
    stop classification after the first database hit; use `--quick`
    to enable this mode.

* **Hit group threshold**: The option `--minimum-hit-groups` will allow
    you to require multiple hit groups (a group of overlapping k-mers that
    share a common minimizer that is found in the hash table) be found
    before declaring a sequence classified,
    which can be especially useful with custom databases when testing
    to see if sequences either do or do not belong to a particular
    genome.

* **Sequence filtering**: Classified or unclassified sequences can be
    sent to a file for later processing, using the `--classified-out`
    and `--unclassified-out` switches, respectively.

* **Output redirection**: Output can be directed using standard shell
    redirection (`|` or `>`), or using the `--output` switch.

* **Compressed input**: Kraken 2 can handle gzip and bzip2 compressed
    files as input by specifying the proper switch of `--gzip-compressed`
    or `--bzip2-compressed`.

* **Input format auto-detection**: If regular files (i.e., not pipes or device files)
    are specified on the command line as input, Kraken 2 will attempt to
    determine the format of your input prior to classification.
    You can disable this by explicitly specifying
    `--gzip-compressed` or `--bzip2-compressed` as appropriate.
    Note that use of the character device file `/dev/fd/0` to read
    from standard input (aka `stdin`) will **not** allow auto-detection.

* **Paired reads**: Kraken 2 provides an enhancement over Kraken 1 in its
    handling of paired read data.  Rather than needing to concatenate the
    pairs together with an `N` character between the reads, Kraken 2 is
    able to process the mates individually while still recognizing the
    pairing information.  Using the `--paired` option to `kraken2` will
    indicate to `kraken2` that the input files provided are paired read
    data, and data will be read from the pairs of files concurrently.

    Usage of `--paired` also affects the `--classified-out` and
    `--unclassified-out` options; users should provide a `#` character
    in the filenames provided to those options, which will be replaced
    by `kraken2` with "`_1`" and "`_2`" with mates spread across the two
    files appropriately.  For example:

        kraken2 --paired --classified-out cseqs#.fq seqs_1.fq seqs_2.fq

    will put the first reads from classified pairs in `cseqs_1.fq`, and
    the second reads from those pairs in `cseqs_2.fq`.

To get a full list of options, use `kraken2 --help`.


Output Formats
==============

Standard Kraken Output Format
-----------------------------

Each sequence (or sequence pair, in the case of paired reads) classified
by Kraken 2 results in a single line of output.  Kraken 2's output lines
contain five tab-delimited fields; from left to right, they are:

1. "C"/"U": a one letter code indicating that the sequence was either
   classified or unclassified.
2. The sequence ID, obtained from the FASTA/FASTQ header.
3. The taxonomy ID Kraken 2 used to label the sequence; this is 0 if
   the sequence is unclassified.
4. The length of the sequence in bp.  In the case of paired read data,
   this will be a string containing the lengths of the two sequences in
   bp, separated by a pipe character, e.g. "98|94".
5. A space-delimited list indicating the LCA mapping of each $k$-mer in
   the sequence(s).  For example, "562:13 561:4 A:31 0:1 562:3" would
   indicate that:
     - the first 13 $k$-mers mapped to taxonomy ID #562
     - the next 4 $k$-mers mapped to taxonomy ID #561
     - the next 31 $k$-mers contained an ambiguous nucleotide
     - the next $k$-mer was not in the database
     - the last 3 $k$-mers mapped to taxonomy ID #562

   Note that paired read data will contain a "`|:|`" token in this list
   to indicate the end of one read and the beginning of another.

   When Kraken 2 is run against a protein database (see [Translated Search]),
   the LCA hitlist will contain the results of querying all six frames of
   each sequence.  Reading frame data is separated by a "`-:-`" token.

Kraken 1 offered a `kraken-translate` and `kraken-report` script to change
the output into different formats.  Through the use of `kraken2 --use-names`,
Kraken 2 will replace the taxonomy ID column with the scientific name and
the taxonomy ID in parenthesis (e.g., "Bacteria (taxid 2)" instead of "2"),
yielding similar functionality to Kraken 1's `kraken-translate` script.
The sample report functionality now exists as part of the `kraken2` script,
with the use of the `--report` option; the sample report formats are
described below.

Sample Report Output Format
---------------------------

Like Kraken 1, Kraken 2 offers two formats of sample-wide results.
Kraken 2's standard sample report format is tab-delimited with one
line per taxon.  The fields of the output, from left-to-right, are
as follows:

1. Percentage of fragments covered by the clade rooted at this taxon
2. Number of fragments covered by the clade rooted at this taxon
3. Number of fragments assigned directly to this taxon
4. A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom,
   (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies.
   Taxa that are not at any of these 10 ranks have a rank code that is
   formed by using the rank code of the closest ancestor rank with
   a number indicating the distance from that rank.  E.g., "G2" is a
   rank code indicating a taxon is between genus and species and the
   grandparent taxon is at the genus rank.
5. NCBI taxonomic ID number
6. Indented scientific name

The scientific names are indented using space, according to the tree
structure specified by the taxonomy.

By default, taxa with no reads assigned to (or under) them will not have
any output produced. However, if you wish to have all taxa displayed, you
can use the `--report-zero-counts` switch to do so. This can be useful if
you are looking to do further downstream analysis of the reports, and want
to compare samples. Sorting by the taxonomy ID (using `sort -k5,5n`) can
provide a consistent line ordering between reports.

In addition, we also provide the option `--use-mpa-style` that can be used
in conjunction with `--report`.  This option provides output in a format
similar to MetaPhlAn's output.  The output with this option provides one
taxon per line, with a lowercase version of the rank codes in Kraken 2's
standard sample report format (except for 'U' and 'R'), two underscores,
and the scientific name of the taxon (e.g., "d__Viruses").  The full
taxonomy of each taxon (at the eight ranks considered) is given, with each
rank's name separated by a pipe character (e.g., "d__Viruses|o_Caudovirales").
Following this version of the taxon's scientific name is a tab and the
number of fragments assigned to the clade rooted at that taxon.


Translated Search
=================

Kraken 2 allows users to perform a six-frame translated search, similar
to the well-known BLASTX program.  To do this, Kraken 2 uses a reduced
15 amino acid alphabet and stores amino acid minimizers in its database.
LCA results from all 6 frames are combined to yield a set of LCA hits,
which is then resolved in the same manner as in Kraken's normal operation.

To build a protein database, the `--protein` option should be given to
`kraken2-build` (either along with `--standard`, or with all steps if
building a custom database).


Custom Databases
================

We realize the standard database may not suit everyone's needs.  Kraken 2
also allows creation of customized databases.

To build a custom database:

 1. Install a taxonomy.  Usually, you will just use the NCBI taxonomy,
    which you can easily download using:

        kraken2-build --download-taxonomy --db $DBNAME

    This will download the accession number to taxon maps, as well as the
    taxonomic name and tree information from NCBI.  These files can
    be found in `$DBNAME/taxonomy/` .  If you need to modify the taxonomy,
    edits can be made to the `names.dmp` and `nodes.dmp` files in this
    directory; you may also need to modify the `*.accession2taxid` files
    appropriately.

    Some of the standard sets of genomic libraries have taxonomic information
    associated with them, and don't need the accession number to taxon maps
    to build the database successfully.  These libraries include all those
    available through the `--download-library` option (see next point), except
    for the `plasmid` and non-redundant databases.  If you are not using
    custom sequences (see the `--add-to-library` option) and are not using
    one of the `plasmid` or non-redundant database libraries, you may want to
    skip downloading of the accession number to taxon maps.  This can be done
    by passing `--skip-maps` to the `kraken2-build --download-taxonomy` command.

 2. Install one or more reference libraries.  Several sets of standard
    genomes/proteins are made easily available through `kraken2-build`:

    - `archaea`: RefSeq complete archaeal genomes/proteins
    - `bacteria`: RefSeq complete bacterial genomes/proteins
    - `plasmid`: RefSeq plasmid nucleotide/protein sequences
    - `viral`: RefSeq complete viral genomes/proteins
    - `human`: GRCh38 human genome/proteins
    - `fungi`: RefSeq complete fungal genomes/proteins
    - `plant`: RefSeq complete plant genomes/proteins
    - `protozoa`: RefSeq complete protozoan genomes/proteins
    - `nr`: NCBI non-redundant protein database
    - `nt`: NCBI non-redundant nucleotide database
    - `UniVec`: NCBI-supplied database of vector, adapter, linker, and
      primer sequences that may be contaminating sequencing projects and/or
      assemblies
    - `UniVec_Core`: A subset of UniVec chosen to minimize false positive
      hits to the vector database

    To download and install any one of these, use the `--download-library`
    switch, e.g.:

        kraken2-build --download-library bacteria --db $DBNAME

    Multiple libraries can be downloaded into a database prior to building
    by issuing multiple `kraken2-build --download-library` commands, e.g.:

        kraken2-build --download-library archaea --db $DBNAME
        kraken2-build --download-library viral --db $DBNAME

    The above commands would prepare a database that would contain archaeal
    and viral genomes; the `--build` option (see below) will still need to
    be used after downloading these libraries to actually build the database,
    however.

    (Note that downloading `nr` requires use of the `--protein`
      option, and that `UniVec` and `UniVec_Core` are incompatible with
      the `--protein` option.)

    Other genomes can also be added, but such genomes must meet certain
    requirements:
    - Sequences must be in a FASTA file (multi-FASTA is allowed)
    - Each sequence's ID (the string between the `>` and the first
      whitespace character on the header line) must contain either
      an NCBI accession number to allow Kraken 2 to lookup the correct taxa,
      or an explicit assignment of the taxonomy ID using `kraken:taxid`
      (see below).

    Sequences not downloaded from NCBI may need their taxonomy information
    assigned explicitly.  This can be done using the string `kraken:taxid|XXX`
    in the sequence ID, with `XXX` replaced by the desired taxon ID.  For
    example, to put a known adapter sequence in taxon 32630 ("synthetic
    construct"), you could use the following:

        >sequence16|kraken:taxid|32630  Adapter sequence
        CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA

    The `kraken:taxid` string must begin the sequence ID or be immediately
    preceded by a pipe character (`|`).  Explicit assignment of taxonomy IDs
    in this manner will override the accession number mapping provided by NCBI.

    If your genomes meet the requirements above, then you can add each
    sequence to your database's genomic library using the `--add-to-library`
    switch, e.g.:

        kraken2-build --add-to-library chr1.fa --db $DBNAME
        kraken2-build --add-to-library chr2.fa --db $DBNAME

    Note that if you have a list of files to add, you can do something like
    this in `bash`:

        for file in chr*.fa
        do
            kraken2-build --add-to-library $file --db $DBNAME
        done

    Or even add all `*.fa` files found in the directory `genomes`:

    `find genomes/ -name '*.fa' -print0 | xargs -0 -I{} -n1 kraken2-build --add-to-library {} --db $DBNAME`

    (You may also find the `-P` option to `xargs` useful to add many files in
    parallel if you have multiple processors.)

3. Once your library is finalized, you need to build the database.  This
   can be done with the command:

        kraken2-build --build --db $DBNAME

   The `--threads` option is also helpful here to reduce build time.

   By default, the values of $k$ and $\ell$ are 35 and 31, respectively (or
   15 and 12 for protein databases).  These values can be explicitly set
   with the `--kmer-len` and `--minimizer-len` options, however.  Note that
   the minimizer length must be no more than 31 for nucleotide databases,
   and 15 for protein databases.  Additionally, the minimizer length $\ell$
   must be no more than the $k$-mer length.  There is no upper bound on
   the value of $k$, but sequences less than $k$ bp in length cannot be
   classified.

   Kraken 2 also utilizes a simple spaced seed approach to increase
   accuracy.  A number $s$ < $\ell$/4 can be chosen, and $s$ positions
   in the minimizer will be masked out during all comparisons.
   Masked positions are chosen to alternate from the second-to-last
   position in the minimizer; e.g., $s$ = 5 and $\ell$ = 31 will result
   in masking out the 0 positions shown here:

       111 1111 1111 1111 1111 1101 0101 0101

   By default, $s$ = 7 for nucleotide databases, and $s$ = 0 for
   protein databases.  This can be changed using the `--minimizer-spaces`
   option along with the `--build` task of `kraken2-build`.

A full list of options for `kraken2-build` can be obtained using
`kraken2-build --help`.

After building a database, if you want to reduce the disk usage of
the database, you can use the `--clean` option for `kraken2-build`
to remove intermediate files from the database directory.

Masking of Low-complexity Sequences
===================================

Low-complexity sequences, e.g. "ACACACACACACACACACACACACAC", are known
to occur in many different organisms and are typically less informative
for use in alignments; the BLAST programs often mask these sequences by
default.  Using this masking can help prevent false positives in Kraken 2's
results, and so we have added this functionality as a default option to
Kraken 2's library download/addition process.

Kraken 2 uses two programs to perform low-complexity sequence masking,
both available from NCBI: `dustmasker`, for nucleotide sequences, and
`segmasker`, for amino acid sequences.  These programs are available
as part of the NCBI BLAST+ suite.  If these programs are not installed
on the local system and in the user's PATH when trying to use
`kraken2-build`, the database build will fail.  Users who do not wish to
install these programs can use the `--no-masking` option to `kraken2-build`
in conjunction with any of the `--download-library`, `--add-to-library`, or
`--standard` options; use of the `--no-masking` option will skip masking of
low-complexity sequences during the build of the Kraken 2 database.

Special Databases
=================

To support some common use cases, we provide the ability to build Kraken 2
databases using data from various external databases.  These external
databases may not follow the NCBI taxonomy, and so we've provided
mechanisms to automatically create a taxonomy that will work with Kraken 2
(although such taxonomies may not be identical to NCBI's).

To build one of these "special" Kraken 2 databases, use the following command:

    kraken2-build --db $DBNAME --special TYPE

where the `TYPE` string is one of the database names listed below.

At present, the "special" Kraken 2 database support we provide is limited
to pre-packaged solutions for some public 16S sequence databases, but this may
grow in the future.

16S Databases
-------------

For targeted 16S sequencing projects, a normal Kraken 2 database using whole
genome data may use more resources than necessary.  A Kraken 2 database created
from a well-curated genomic library of just 16S data can provide both a more
efficient solution as well as a more accurate set of predictions for such
projects.  We provide support for building Kraken 2 databases from three
publicly available 16S databases:

* [Greengenes] (Kraken 2 database name: `greengenes`), using all available 16S data.
* [RDP] (Kraken 2 database name: `rdp`), using the bacterial and archaeal 16S data.
* [SILVA] (Kraken 2 database name: `silva`), using the Small subunit NR99 sequence set.

Note that these databases may have licensing restrictions regarding their data,
and it is your responsibility to ensure you are in compliance with those
restrictions; please visit the databases' websites for further details.  The
`kraken2-build` script only uses publicly available URLs to download data and
then converts that data into a form compatible for use with Kraken 2.

Furthermore, if you use one of these databases in your research, please
visit the corresponding database's website to determine the appropriate and
up-to-date citation.

[Greengenes]:   http://greengenes.lbl.gov/
[RDP]:          http://rdp.cme.msu.edu/
[SILVA]:        http://www.arb-silva.de/

Confidence Scoring
==================

At present, we have not yet developed a confidence score with a
probabilistic interpretation for Kraken 2.  However, we have developed a
simple scoring scheme that has yielded good results for us, and we've
made that available in Kraken 2 through use of the `--confidence` option
to `kraken2`.  The approach we use allows a user to specify a threshold
score in the [0,1] interval; the classifier then will adjust labels up
the tree until the label's score (described below) meets or exceeds that
threshold.  If a label at the root of the taxonomic tree would not have
a score exceeding the threshold, the sequence is called unclassified by
Kraken 2 when this threshold is applied.

A sequence label's score is a fraction $C$/$Q$, where $C$ is the number of
$k$-mers mapped to LCA values in the clade rooted at the label, and $Q$ is the
number of $k$-mers in the sequence that lack an ambiguous nucleotide (i.e.,
they were queried against the database).  Consider the example of the
LCA mappings in Kraken 2's output given earlier:

"562:13 561:4 A:31 0:1 562:3" would indicate that:

* the first 13 $k$-mers mapped to taxonomy ID #562
* the next 4 $k$-mers mapped to taxonomy ID #561
* the next 31 $k$-mers contained an ambiguous nucleotide
* the next $k$-mer was not in the database
* the last 3 $k$-mers mapped to taxonomy ID #562

In this case, ID #561 is the parent node of #562.  Here, a label of #562
for this sequence would have a score of $C$/$Q$ = (13+3)/(13+4+1+3) = 16/21.
A label of #561 would have a score of $C$/$Q$ = (13+4+3)/(13+4+1+3) = 20/21.
If a user specified a `--confidence` threshold over 16/21, the classifier
would adjust the original label from #562 to #561; if the threshold was
greater than 20/21, the sequence would become unclassified.

Inspecting a Kraken 2 Database's Contents
=========================================

The `kraken2-inspect` script allows users to gain information about the content
of a Kraken 2 database.  The output format of `kraken2-inspect`
is identical to the reports generated with the `--report` option to `kraken2`.
Instead of reporting how many reads in input data classified to a given taxon
or clade, as `kraken2`'s `--report` option would, the `kraken2-inspect` script
will report the number of minimizers in the database that are mapped to the
various taxa/clades.  For example, the first five lines of `kraken2-inspect`'s
output on an example database might look like this:

    $ kraken2-inspect --db EXAMPLE_DB | head -5
    100.00% 1770368409      1581179 R       1       root
     96.50% 1708407622      58003   R1      131567    cellular organisms
     91.28% 1615910070      985309  D       2           Bacteria
     43.89% 777062062       1312736 P       1224          Proteobacteria
     18.62% 329590216       555667  C       1236            Gammaproteobacteria

This output indicates that 555667 of the minimizers in the database map
directly to the Gammaproteobacteria class (taxid #1236), and 329590216 (18.62%)
of the database's minimizers map to a taxon in the clade rooted at
Gammaproteobacteria.  For more information on `kraken2-inspect`'s options,
use its `--help` option.


Distinct minimizer count information
====================================

The [KrakenUniq] project extended Kraken 1 by, among other things, reporting
an estimate of the number of distinct k-mers associated with each taxon in the
input sequencing data.  This allows users to better determine if Kraken's
classifications are due to reads distributed throughout a reference genome,
or due to only a small segment of a reference genome (and therefore likely
false positive).

Thanks to the generosity of KrakenUniq's developer Florian Breitwieser in
allowing parts of the KrakenUniq source code to be licensed under Kraken 2's
MIT license, this distinct counting estimation is now available in Kraken 2.
Development work by Martin Steinegger and Ben Langmead helped bring this
functionality to Kraken 2.

At present, this functionality is an optional *experimental feature* -- meaning
that we may later alter it in a way that is not backwards compatible with
previous versions of the feature.

To use this functionality, simply run the `kraken2` script with the additional
`--report-minimizer-data` flag along with `--report`, e.g.:

    kraken2 --db $DBNAME --report k2_report.txt --report-minimizer-data \
        --output k2_output.txt sequence_data.fq

This will put the standard Kraken 2 output (formatted as described in
[Standard Kraken Output Format]) in `k2_output.txt` and the report information
in `k2_report.txt`.  Within the report file, two additional columns will be
present, e.g.:

**normal report format**:

    36.40	182	182	S2	211044	                      Influenza A virus (A/Puerto Rico/8/1934(H1N1))

**modified report format**:

    36.40	182	182	1688	18	S2	211044	                      Influenza A virus (A/Puerto Rico/8/1934(H1N1))

In this modified report format, the two new columns are the fourth and fifth,
respectively representing the number of minimizers found to be associated with
a taxon in the read sequences (1688), and the estimate of the number of distinct
minimizers associated with a taxon in the read sequence data (18).  This would
indicate that although 182 reads were classified as belonging to H1N1 influenza,
only 18 distinct minimizers led to those 182 classifications.

The format with the `--report-minimizer-data` flag, then, is similar to that
described in [Sample Report Output Format], but slightly different.  The fields
in this new format, from left-to-right, are:

1. Percentage of fragments covered by the clade rooted at this taxon
2. Number of fragments covered by the clade rooted at this taxon
3. Number of fragments assigned directly to this taxon
4. Number of minimizers in read data associated with this taxon (**new**)
5. An estimate of the number of distinct minimizers in read data associated
   with this taxon (**new**)
6. A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom,
   (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies.
   Taxa that are not at any of these 10 ranks have a rank code that is
   formed by using the rank code of the closest ancestor rank with
   a number indicating the distance from that rank.  E.g., "G2" is a
   rank code indicating a taxon is between genus and species and the
   grandparent taxon is at the genus rank.
7. NCBI taxonomic ID number
8. Indented scientific name

We decided to make this an optional feature so as not to break existing
software that processes Kraken 2's standard report format.  However, this
new format can be converted to the standard report format with the command:

    cut -f1-3,6-8 k2_new_report.txt > k2_std_report.txt

As noted above, this is an *experimental feature*.  We intend to continue
development on this feature, and may change the new format and/or its
information if we determine it to be necessary.

For background on the data structures used in this feature and their
interaction with Kraken, please read the [KrakenUniq paper], and please
cite that paper if you use this functionality as part of your work.

[KrakenUniq]: https://github.com/fbreitwieser/krakenuniq
[KrakenUniq paper]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1568-0


Kraken 2 Environment Variables
==============================

The `kraken2` and `kraken2-inspect` scripts supports the use of some
environment variables to help in reducing command line lengths:

* **`KRAKEN2_NUM_THREADS`**: if the
    `--threads` option is not supplied to `kraken2`, then the value of this
    variable (if it is set) will be used as the number of threads to run
    `kraken2`.  (This variable does not affect `kraken2-inspect`.)

* **`KRAKEN2_DB_PATH`**: much like the `PATH` variable is used for executables
    by your shell, `KRAKEN2_DB_PATH` is a colon-separated list of directories
    that will be searched for the database you name if the named database
    does not have a slash (`/`) character.  By default, Kraken 2 assumes the
    value of this variable is "`.`" (i.e., the current working directory).
    This variable can be used to create one (or more) central repositories
    of Kraken databases in a multi-user system.  Example usage in bash:

        export KRAKEN2_DB_PATH="/home/user/my_kraken2_dbs:/data/kraken2_dbs:"

    This will cause three directories to be searched, in this order:

    1) `/home/user/my_kraken2_dbs`
    2) `/data/kraken2_dbs`
    3) the current working directory (caused by the empty string as
         the third colon-separated field in the `KRAKEN2_DB_PATH` string)

    The search for a database will stop when a name match is found; if
    two directories in the `KRAKEN2_DB_PATH` have databases with the same
    name, the directory of the two that is searched first will have its
    database selected.

    If the above variable and value are used, and the databases
    `/data/kraken2_dbs/mainDB` and `./mainDB` are present, then

        kraken2 --db mainDB sequences.fa

    will classify `sequences.fa` using `/data/kraken_dbs/mainDB`; if instead
    you wanted to use the `mainDB` present in the current directory,
    you would need to specify a directory path to that database in order
    to circumvent searching, e.g.:

        kraken2 --db ./mainDB sequences.fa

    Note that the `KRAKEN2_DB_PATH` directory list can be skipped by the use
    of any absolute (beginning with `/`) or relative pathname (including
    at least one `/`) as the database name.

* **`KRAKEN2_DEFAULT_DB`**: if no database is supplied with the `--db` option,
    the database named in this variable will be used instead.  Using this
    variable, you can avoid using `--db` if you only have a single database
    that you usually use, e.g. in bash:

        export KRAKEN2_DEFAULT_DB="/home/user/kraken2db"
        kraken2 sequences.fa > kraken2.output

    This will classify `sequences.fa` using the `/home/user/kraken2db`
    database.

    Note that the value of `KRAKEN2_DEFAULT_DB` will also be interpreted in
    the context of the value of `KRAKEN2_DB_PATH` if you don't set
    `KRAKEN2_DEFAULT_DB` to an absolute or relative pathname.  Given the earlier
    example in this section, the following:

        export KRAKEN2_DEFAULT_DB="mainDB"
        kraken2 sequences.fa

    will use `/data/kraken_dbs/mainDB` to classify `sequences.fa`.
