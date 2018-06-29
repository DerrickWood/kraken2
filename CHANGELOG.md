# Changelog

## Unreleased

### Added
- Disk usage info to kraken2-build --clean

### Changed
- Move to /usr/bin/env for perl scripts
- Add DB loading message to keep people from killing processes early
- Add flag files for resuming download of nucleotide accession map data

### Fixed
- Allow d/l of protozoa library w/ kraken2-build script
- Filenames for SILVA database taxonomy info

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
