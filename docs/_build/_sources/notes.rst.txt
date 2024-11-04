=====
Notes
=====

* The `genomics` DataFrame accessor

  * Will only work if the entire DataFrame consists of GenotypeArray Series
  * Requires that all variant IDs are unique.  Variants get a random unique (UUID4) ID if one is not specified.

* The Series (or DataFrame column) name should not be confused with the variant ID.  There is no reason to assume they match.
