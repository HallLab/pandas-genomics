===============
Release History
===============

v0.11.0 (2021-07-21)
--------------------

Add a *Region* scalar with support for loading from BED files and filtering

v0.10.1 (2021-07-16)
--------------------

Fixed an issue reporting the size of GenotypeDtype in memory

v0.10.0 (2021-07-16)
--------------------

Change genotype scores to be a unit8 (255=missing) instead of float64, to save ~70% of memory usage

v0.9.1 (2021-07-13)
-------------------

Fix a side-effect of to_plink that swapped allele idxs which could affect edge encoding

v0.9.0 (2021-07-12)
-------------------

Remove fancy new pandas dtypes from edge encoding (and rename weighted encoding to edge encoding).

v0.8.0 (2021-07-09)
-------------------

Improve randomness handling in BAMS and fix calculations to ensure they work as expected (especially quantitative)

v0.7.0 (2021-07-06)
-------------------

Return float values for encoding instead of newer (and less compatible) pandas types

v0.6.1 (2021-06-26)
-------------------

Specify only alternate allele frequencies when generating random genotypes, to avoid float rounding problems

v0.6.0 (2021-06-25)
-------------------

Enhancements
^^^^^^^^^^^^

* The *genomics* DataFrame Accessor no longer requires that all columns in the DataFrame are backed by a GenotypeArray

v0.5.2 (2021-06-24)
-------------------

* Update numpy version requirement

v0.5.1 (2021-06-23)
-------------------

* Update dependencies to fix doc building

v0.5.0 (2021-06-23)
-------------------

Enhancements
^^^^^^^^^^^^
* Add separate Series and DataFrame Accessors
* Speed up Plink file import
* Add prototype BiAllelic Model Simulator
* Add random genotype simulator
* Add Plink Output

v0.4.0 (2021-03-30)
-------------------

Enhancements
^^^^^^^^^^^^
* Updates to the backing structure of the GenotypeArray
* Added Variant quality scores
* Added Genotype quality scores

v0.3.0 (2020-10-23)
-------------------

Enhancements
^^^^^^^^^^^^
* Four encoding methods added to `GenotypeArray`: `encode_additive`, `encode_recessive`, `encode_dominant`,
  and `encode_codominant`
* A `genotype` Series accessor has been added.  Thus far it includes the encoding methods.
* Variants now track which allele is considered the reference allele

Tests
^^^^^
* I/O Tests are added for Plink and VCF.  Note that VCF requires htslib, so tests are skipped on Windows
* Encoding tests added

Docs
^^^^
Improved.

v0.2.0 (2020-10-09)
-------------------

Tests
^^^^^
Pandas ExtensionArray tests are passing.
Note: a few can't pass until it's possible to define dtype as scalar
(See `Pandas #33825  <https://github.com/pandas-dev/pandas/issues/33825>`_)

Docs
^^^^
Some initial documentation added

v0.1.0 (2020-09-18)
-------------------

Initial Release