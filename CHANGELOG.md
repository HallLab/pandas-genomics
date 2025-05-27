# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.1] – 2025-05-27

### Added
- Full compatibility with `pandas.tests.extension.base` test suite.
- `__setitem__`: added support for `Genotype`, `GenotypeArray`, and `pd.Series` (with `GenotypeArray` values).

### Changed
- `__getitem__`: uses `pandas.api.indexers.check_array_indexer` and raises clear `ValueError` for invalid keys like strings.
- `_from_sequence`: improved validation logic for variant compatibility and handling of scalar values.
- `factorize`: now properly detects unique values using `allele_idxs` and handles `na_value`.

### Fixed
- `insert()`: properly rejects non-Genotype scalars and raises informative errors.
- `__eq__`: avoids invalid comparisons by checking variant compatibility.
- Compatibility adjustments for `pandas 2.x`, `Python 3.11+`.

### Deprecated
- None.

### Removed
- None.

### Security
- None.
