# Changelog

## [2.1.6] - 2025-07-08

### Changed
- Updated GTDB version from 220 to 226

### Added
- Added BLAST files to FASTA convertor, this currently only works for nucleotides
- Added support to `k2` for downloading and processing `core_nt`, `nt`, `env_nt` and `viral_nt` BLAST databases
- Added `--blast-volumes` option for allowing users the ability to download specific volumes
- Added `--gtdb-server` option to `k2 build` allowing users to select the appropriate GTDB server

### Fixed
- Fixed an issue causing undefined behavior in `k2` due to sharing of global variables between processes
- Fixed an issue causing `k2` to hang when processing compressed reads
- Fixed an issue causing "classify daemon" to behave incorrectly when processing compressed reads
- Fixed an issue causing `classify` to crash when using the `--memory-mapping` option
- Fixed an issue preventing `k2` from fetching assemblies from both `refseq` and `genbank` when 
  `--assembly-source` is set to `all` (@beantkapoor786: [PR 968](https://github.com/DerrickWood/kraken2/pull/968))

## [2.1.5] - 2025-04-18

### Added
- Added experimental support for running the classifier as a daemon. The classifier
  daemon will persist any database that it has loaded allowing users to save time
  when running multiple samples on the same large database. This feature is currently
  only available with the `k2` wrapper by specifying the `--use-daemon` option when
  running `k2 classify`.

  Example:
  `k2 classify --db pluspfp --use-daemon --threads 12 sample.fa --output sample.out`
  On first invocation the pluspfp database will be loaded into memory. On subsequent runs of
  this command the database load step will be skipped. If a new database is provided that
  database will persisted in memory as well.

- Added `--stop-daemon` option to `k2 clean`. This option will stop a running classifier
  daemon.
- `--resume` option has been added to `k2 download-library`. `k2 download-library` has
  been updated so that files downloaded from NCBI will be saved with a `.tmp` extension. This
  extension will only be removed after `k2` has verified that the MD5 sum of the downloaded
  file matches that of the server. The `--resume` option builds off of that foundation by
  skipping files that do not have the `.tmp` extension. If a user has a library that was
  downloaded with a version of kraken 2 prior to version 2.1.5, the resume option will fail if
  partially downloaded files exist. Run `k2 download-library` without the `--resume` option
  first to ensure that all files are fully downloaded after which resume can be used successfully.
  When using 2.1.5 to download a new library, `--resume` should work without issue on a failed
  download.

### Changed
- Made general improvements to `k2` wrapper
- Removed many threads option from `dump_table`
- `--assembly-level` option now works with library collections

### Fixed
- Fixed memory regressions in Kraken 2
- Fixed an issue causing `classify` to incorrectly format unclassified sequences
- Fixed an issue causing zero-counts to sometimes have missing entries
- Fixed an issue causing Kraken 2 to sometimes seg-fault when inspecting files
- Fixed an issue causing `k2mask` to read only from stdin
- Fixed an issue causing `classify` to not process paired reads

## [2.1.4] - 2025-02-17

### Added
- Added support for building GTDB databases
- Added memory-mapping support to dump<sub>table</sub> (kraken2-inspect)
- Added support for many threads to dump<sub>table</sub> (kraken2-inspect)
- Added support for downloading genome accessions by project ID, taxon ID or accession number to `k2`.
  This feature makes use of the NCBI API <https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/>

### Changes
- Replaced FASTA/Q parser with kseq resulting in significantly faster runtimes when estimating database capacity
- Sped up decompression in `k2` by using `zlib` instead of `gzip` making `k2` competitive with native gzip binaries
- Changed the way accessions are downloaded when using NCBI API so that downloads are more reliable
- Significant improvements to stability and speed of `k2`
- Updated SILVA from 138.1 to 138.2

## [2.1.3] - 2023-06-06

### Added
- k2 python wrapper script

### Changed
- Low complexity masking for nucleotide sequences now performed by our own
  multithreaded code (k2mask) instead of the dustmasker program

### Fixed
- Bug that caused infinite loop when sequences had taxid 0 (thanks to R. Charles)
- Silenced new warning in GCC 11
- Allow NCBI downloads from both ftp:// and https:// paths (thanks to M. Machado)

## [2.1.2] - 2021-05-10

### Changed
- New retry behavior for kraken2-build --download-library --use-ftp

### Fixed
- Bug that caused classification to not work when k-mer length was equal
  to minimizer length
- Missing source file in CMakeLists.txt

## [2.1.1] - 2020-11-08

### Fixed
- Compilation error with GCC 10, omission of cstdint header
- Removed --skip-maps from standard install due to addition of plasmids,
  which require acc/taxid maps

## [2.1.0] - 2020-10-13

### Added
- Small viral reference set and read simulator for future testing
- Build options to aid testing
- Integration of HLL to estimate distinct minimizer counts

### Changed
- Build code now creates databases with deterministic MD5 sums by default;
  --fast-build option to kraken2-build introduced to access old behavior
- Added plasmid library to standard installation set
- Updated SILVA to release 138.1

### Fixed
- Modified build code to prevent insertion of minimizers with ambiguous bases
- Bug where hit list output in quick mode used internal taxid (not external)
- No more attempts to download "na" paths from FTP site
- Runaway memory usage bug with unpaired classification (thanks to D. Cameron)
- Modified estimation code to better handle small reference libraries

## [2.0.9] - 2020-04-07 (beta)

### Added
- Expose --load-factor setting to kraken2-build
- New --minimum-hit-groups option to kraken2

### Changed
- Require 2 hit groups (set of overlapping k-mers w/ same minimizer) to
  make classification by default
- Allow build options to pass through to subsequent invocations (e.g.,
  k-mer length for 16S DBs)
- Removed env options for library downloads (no longer available from
  same NCBI location)
- Updated SILVA to release 138

### Fixed
- Removed mention of --fastq-input from Manual
- Made PE read identifier suffix trimming more restrictive (only on /1 and /2)
- Bug where some reads would be classified with taxid 0
- Bug that didn't allow kraken2-inspect to work with large databases

## [2.0.8] - 2019-04-25 (beta)

### Added
- FTP downloading option for taxonomy/libraries (--use-ftp for kraken2-build)
- Option to skip downloading taxonomy maps

### Changed
- Added lookup table to speed up parsing in MinimizerScanner class
- Default parameters for minimizer lengths and spaces changed (spaces=7 for
  nucleotide search, length=12 for translated search)

### Fixed
- Linked space expansion value for proteins to constant used by MinimizerScanner
- Reporting of taxids in classified-out sequence files
- Confidence scoring bug associated with failure to leave some sequences
  unclassified
- Reverse complement shifting bug, code made backwards-compatible with
  existing databases (newly created DBs will have fix)
- NCBI taxonomy download error due to removal of EST/GSS files

## [2.0.7] - 2018-08-11 (beta)

### Added
- Disk usage info to kraken2-build --clean
- Memory allocation error message for hash table
- Option for --max-db-size and hash downsampling
- Multithreading to kraken2-inspect

### Changed
- Move to /usr/bin/env for perl scripts
- Add DB loading message to keep people from killing processes early
- Add flag files for resuming download of nucleotide accession map data
- Converted lookup_accession_numbers script into C++ program w/ memory mapping
- Clarified in manual that one or more libraries allowed for custom DBs
- Silenced progress messages in C++ programs for non-TTY stderr
- Taxonomy downloads switched to rsync from wget (ftp)
- Removed '%' from reports

### Fixed
- Allow d/l of protozoa library w/ kraken2-build script
- Filenames for SILVA database taxonomy info
- Typo in manual for output format example
- Corrected default space count in manual
- Removed obvious race condition in --add-to-library functionality
- Corrected behavior of --classified-out and --unclassified-out (no longer
  forcing .fq/.fa file extensions, respecting '#' in paired mode)
- Usage message in kraken2-inspect
- Taxonomy creation for 16S databases

## [2.0.6] - 2018-06-13 (beta)

### Changed
- New DB summary info printed out w/ inspect script + --skip-counts option

### Fixed
- Now stripping carriage returns and other trailing whitespace from sequence
  data
- Treating l-mers immediately following ambiguous characters as ambiguous
  until a full k-mer is processed
- Bug in expansion of spaced seed masks that left spaces at end

## [2.0.5] - 2018-05-21 (beta)

### Added
- New kraken2-inspect script to report minimizer counts per taxon

### Changed
- Kraken 2X build now adds terminators to all reference sequences

## [2.0.4] - 2018-05-06 (beta)

### Fixed
- Improved portability to older g++ by removing initialization of
  variable-length string.

## [2.0.3] - 2018-02-12 (alpha)

### Added
- Reporting options to kraken2 script (like Kraken 1's kraken-report and
  kraken-mpa-report)

### Changed
- Made loading to RAM default option, added --memory-mapping option to kraken2

## [2.0.2] - 2018-02-04 (alpha)

### Added
- Low base quality masking option

### Changed
- Moved low-complexity masking to library download/addition, out of build
  process
- Made no masking default for human genome in standard installation

## [2.0.1] - 2018-01-01 (alpha)

### Added
- Low-complexity sequence masking as a default
- UniVec/UniVec_Core databases to supported downloads
- UniVec_Core & human in standard Kraken 2 DB
- 16S DB support (Greengenes, Silva, RDP)
- --use-names flag for kraken2 script
- Priority queue to ensure classifier output order matches input order when
  multi-threading
- Changelog

### Changed
- Reduced amino acid alphabet (requires rebuild of old protein DBs)
- Operating manual

### Fixed
- kraken2 now allows compression & paired processing at same time

## [2.0.0] - 2017-12-04 (alpha, initial release)
