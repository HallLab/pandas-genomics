===============
Release History
===============

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