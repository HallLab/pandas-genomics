.. pandas-genomics documentation master file, created by
   sphinx-quickstart on Mon Sep 21 13:03:41 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for pandas-genomics
=================================
:Version: |version|

Pandas ExtensionArray for working with genomics data

Quickstart
----------

`Variant` objects holds information about a particular variant:

.. code-block:: python

    >>> from pandas_genomics.scalars import Variant
    >>> variant = Variant('12', 112161652, id='rs12462', ref='A', alt=['C', 'T'], score=30)
    >>> print(variant)
    rs12462(n2)[chr=12;pos=112161652;ref=A;alt=C,T]Q30

`Genotype` objects are associated with a particular `Variant`:

.. code-block:: python

    >>> gt = variant.make_genotype("A", "C")
    >>> print(gt)
    A/C

The `GenotypeArray` stores genotypes with an associated variant and has useful methods and properties.
The `genotype` Series accessor allows those methods and properties to be accessed from a Series.

.. code-block:: python

    >>> from pandas_genomics.scalars import Variant
    >>> from pandas_genomics.arrays import GenotypeArray
    >>> variant = Variant('12', 112161652, id='rs12462', ref='A', alt=['C'])
    >>> gt_array = GenotypeArray([variant.make_genotype_from_str(s) for s in ["C/C", "A/C", "A/A"]])

    >>> print(gt_array)
    <GenotypeArray>
    [Genotype(variant=rs12462[chr=12;pos=112161652;ref=A;alt=C])[C/C],
     Genotype(variant=rs12462[chr=12;pos=112161652;ref=A;alt=C])[A/C],
     Genotype(variant=rs12462[chr=12;pos=112161652;ref=A;alt=C])[A/A]]
    Length: 3, dtype: genotype(2n)[12; 112161652; rs12462; A; C]

    >>> print(gt_array.astype(str))
    ['C/C' 'A/C' 'A/A']

    >>> print(gt_array.encode_dominant())
    <IntegerArray>
    [1.0, 1.0, 0.0]
    Length: 3, dtype: float

    >>> print(pd.Series(gt_array).genotype.encode_dominant())
    0    1.0
    1    1.0
    2    0.0
    Name: rs12462_C, dtype: float

    >>> print(pd.Series(gt_array).genotype.variant)
    rs12462[chr=12;pos=112161652;ref=A;alt=C]

    >>> print(pd.Series(gt_array).genotype.gt_scores)
    [nan nan nan]

    >>> import pandas as pd
    >>> print(pd.Series(gt_array).genotype.encode_codominant())
    0    Hom
    1    Het
    2    Ref
    Name: rs12462_C, dtype: category
    Categories (3, object): ['Ref' < 'Het' < 'Hom']

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Hall Lab Homepage <https://www.hall-lab.org>
   Pandas-Genomics Github Repo <https://github.com/HallLab/pandas-genomics>

API Reference
-------------

If you are looking for information on a specific function, class or
method, this part of the documentation is for you.

.. toctree::
   :maxdepth: 3

   api

Additional Notes
----------------

Release History, etc

.. toctree::
   :maxdepth: 2

   notes
   release-history