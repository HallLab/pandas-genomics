<div align="center">
<img src="https://github.com/HallLab/pandas-genomics/raw/master/docs/_static/logo.png" alt="pandas_genomics logo"/>
</div>

<br/>

<div align="center">

<!-- Python version -->
<a href="https://pypi.python.org/pypi/pandas-genomics">
<img src="https://img.shields.io/badge/python-3.7+-blue.svg?style=flat-square" alt="PyPI version"/>
</a>
<!-- PyPi -->
<a href="https://pypi.org/project/pandas-genomics/">
<img src="https://img.shields.io/pypi/v/pandas-genomics.svg?style=flat-square" alt="pypi" />
</a><br>
<!-- Build status -->
<a href="https://github.com/HallLab/pandas-genomics/actions?query=workflow%3ACI">
<img src="https://img.shields.io/github/workflow/status/HallLab/pandas-genomics/CI?style=flat-square" alt="Build Status" />
</a>
<!-- Docs -->
<a href="https://pandas-genomics.readthedocs.io/en/latest/">
<img src="https://img.shields.io/readthedocs/pandas-genomics?style=flat-square" alt="Read the Docs" />
</a>
<!-- Test coverage -->
<a href="https://codecov.io/gh/HallLab/pandas-genomics/">
<img src="https://img.shields.io/codecov/c/gh/HallLab/pandas-genomics.svg?style=flat-square" alt="Coverage Status"/>
</a><br>
<!-- License -->
<a href="https://opensource.org/licenses/BSD-3-Clause">
<img src="https://img.shields.io/pypi/l/pandas-genomics?style=flat-square" alt="license"/>
</a>
<!-- Black -->
<a href="https://github.com/psf/black">
<img src="https://img.shields.io/badge/code%20style-Black-black?style=flat-square" alt="code style: black"/>
</a>
</div>

<br/>

Pandas ExtensionDtypes and ExtensionArray for working with genomics data

Quickstart
----------

`Variant` objects holds information about a particular variant:

```python
from pandas_genomics.scalars import Variant
variant = Variant('12', 112161652, id='rs12462', ref='A', alt=['C', 'T'])
print(variant)
```
    rs12462[chr=12;pos=112161652;ref=A;alt=C,T]
    
Each variant should have a unique ID, and a random ID is generated if one is not specified.

`Genotype` objects are associated with a particular `Variant`:

```python
gt = variant.make_genotype("A", "C")
print(gt)
```
```
A/C
```

The `GenotypeArray` stores genotypes with an associated variant and has useful methods and properties:

```python
from pandas_genomics.scalars import Variant
from pandas_genomics.arrays import GenotypeArray
variant = Variant('12', 112161652, id='rs12462', ref='A', alt=['C'])
gt_array = GenotypeArray([variant.make_genotype_from_str(s) for s in ["C/C", "A/C", "A/A"]])
print(gt_array)
```

```
<GenotypeArray>
[Genotype(variant=rs12462[chr=12;pos=112161652;ref=A;alt=C], allele1=1, allele2=1),
Genotype(variant=rs12462[chr=12;pos=112161652;ref=A;alt=C], allele1=0, allele2=1),
Genotype(variant=rs12462[chr=12;pos=112161652;ref=A;alt=C], allele1=0, allele2=0)]
Length: 3, dtype: genotype[12; 112161652; rs12462; A; C]
```

```python
print(gt_array.astype(str))
```

```
    ['C/C' 'A/C' 'A/A']
```

```python
print(gt_array.encode_dominant())
```

```
    <IntegerArray>
    [1.0, 1.0, 0.0]
    Length: 3, dtype: float
```

There are also `genomics` accessors for Series and DataFrame

```python
import pandas as pd
print(pd.Series(gt_array).genomics.encode_codominant())
```

```
    0    Hom
    1    Het
    2    Ref
    Name: rs12462_C, dtype: category
    Categories (3, object): ['Ref' < 'Het' < 'Hom']
```
