import operator
import re
from typing import Dict, MutableMapping, Any, Optional, List, Union, Tuple

import numpy as np
import pandas as pd
from pandas.core.arrays import ExtensionArray
from pandas.core.dtypes.dtypes import register_extension_dtype, PandasExtensionDtype

from .scalars import Variant, Genotype


@register_extension_dtype
class GenotypeDtype(PandasExtensionDtype):
    """
    An ExtensionDtype for genotype data.

    Parameters
    ----------
    variant: Variant
        The ~Variant associated with the genotype

    Attributes
    ----------
    variant

    Examples
    --------
    v = Variant(variant_id="rs12462", chromosome="12", coordinate=112161652, alleles=["T", "C"])
    >>> GenotypeDtype(v)
    genotype[rs12462; 12; 112161652; T,C]
    """

    # Internal attributes
    # -------------------
    # metadata field names
    _metadata = ("variant",)
    # Regular expression
    _match = re.compile(r"(genotype)\["
                        r"(?P<variant_id>.+); "
                        r"(?P<chromosome>.+); "
                        r"(?P<coordinate>[0-9]+); "
                        r"(?P<alleles>.+)\]")
    kind = 'O'
    type = Genotype
    _record_type = np.dtype([('allele1', '>u8'), ('allele2', '>u8')])

    # ExtensionDtype Properties
    # -------------------------
    @property
    def na_value(self) -> Genotype:
        """
        Return the genotype with variant information but no alleles specified
        """
        return Genotype(variant=self.variant)

    @property
    def name(self) -> str:
        return str(self)

    # init
    # ----
    def __init__(self, variant: Variant):
        self.variant = variant

    # ExtensionDtype Methods
    # -------------------------
    @classmethod
    def construct_array_type(cls) -> type:
        """
        Return the array type associated with this dtype

        Returns
        -------
        type
        """
        return GenotypeArray

    @classmethod
    def construct_from_string(cls, string):
        """
        Construct a GenotypeDtype from a string.

        Parameters
        ----------
        string : str
            The string alias for this GenotypeDtype.
            Should be formatted like `genotype[<variant>; <chromosome>; <coordinate>; <comma separated alleles>]`

        Examples
        --------
        >>> GenotypeDtype.construct_from_string('genotype[rs12462; 12; 112161652; T,C]')
        genotype[rs12462; 12; 112161652; T,C]
        """
        if isinstance(string, str):
            msg = "Cannot construct a 'GenotypeDtype' from '{}'"
            try:
                match = cls._match.match(string)
                if match is not None:
                    d = match.groupdict()
                    variant = Variant(variant_id=d["variant_id"],
                                      chromosome=d["chromosome"],
                                      coordinate=int(d["coordinate"]),
                                      alleles=d["alleles"].split(','))
                    return cls(variant=variant)
                else:
                    raise TypeError(msg.format(string))
            except Exception:
                raise TypeError(msg.format(string))
        else:
            raise TypeError(f"'construct_from_string' expects a string, got {type(string)}>")

    @classmethod
    def from_genotype(cls, genotype: Genotype):
        """
        Construct a GenotypeDtype from a Genotype.

        Parameters
        ----------
        genotype: Genotype

        Examples
        -------
        >>> variant = Variant('12', 112161652, 'rs12462')
        >>> genotype = variant.make_genotype_from_str('C/T')
        >>> GenotypeDtype.from_genotype(genotype)
        genotype[rs12462; 12; 112161652; T,C]
        """
        return cls(genotype.variant)

    @classmethod
    def is_dtype(cls, dtype) -> bool:
        """
        Return a boolean if the passed type can be processed as this dtype
        """
        if isinstance(dtype, cls):
            return True
        elif isinstance(dtype, str):
            if dtype.lower().startswith("genotype["):
                try:
                    if cls.construct_from_string(dtype) is not None:
                        return True
                    else:
                        return False
                except (ValueError, TypeError):
                    return False
            else:
                return False
        return super().is_dtype(dtype)

    # Other methods
    # -------------

    def __str__(self):
        return f"genotype[" \
               f"{self.variant.variant_id}; " \
               f"{self.variant.chromosome}; " \
               f"{self.variant.coordinate}; " \
               f"{','.join(self.variant.alleles)}]"

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        if isinstance(other, str):
            return other == self.name

        return (
                isinstance(other, GenotypeDtype)
                and self.variant == other.variant
        )

    # Pickle compatibility
    # --------------------
    def __getstate__(self) -> Dict[str, Any]:
        # pickle support; we don't want to pickle the cache
        return {k: getattr(self, k, None) for k in self._metadata}

    def __setstate__(self, state: MutableMapping[str, Any]) -> None:
        self._variant = state.pop("variant")

    # Other internal methods
    # ----------------------

    def unpack_genotype(self, genotype: Optional[Genotype]):
        if genotype is None:
            return self.unpack_genotype(self.na_value)
        elif not isinstance(genotype, Genotype):
            raise ValueError(f"The fill value must be None or a Genotype object, not {type(genotype)}")
        else:
            return np.array((genotype.allele1, genotype.allele2), dtype=self._record_type)


class GenotypeArray(ExtensionArray):
    """
    Holder for genotypes

    Variant information is stored as part of the type, and the genotype is stored as a pair of integer arrays

    Parameters
    ----------
    values : list-like
        The values of the genotypes.
    dtype : GenotypeDtype
        The specific parametized type.  Optional (if possible to infer from `values`)

    Attributes
    ----------
    dtype: GenotypeDtype
        The specific parametized type
    data: np.dtype([('allele1', '>u8'), ('allele2', '>u8')])
        The genotype values encoded as indices into the allele list of the dtype

    Methods
    -------

    Examples
    --------

    """

    def __init__(self, values: Union[List[Genotype], 'GenotypeArray', np.ndarray],
                 dtype: Optional[GenotypeDtype] = None, copy: bool = False):
        """Initialize assuming values is a GenotypeArray or a numpy array with the correct underlying shape"""
        # If the dtype is passed, ensure it is the correct type
        if GenotypeDtype.is_dtype(dtype):
            self._dtype = dtype
        elif dtype is None:
            self._dtype = None
        else:
            raise ValueError(f"The passed dtype '{dtype}' is not a GenotypeDtype")

        # Load the values
        # ---------------
        if isinstance(values, np.ndarray) and (values.dtype == GenotypeDtype._record_type):
            # Stored data format
            self._data = values

        elif self.is_genotype_array(values):
            # values is a GenotypeArray, simply check the dtype and return
            if self.dtype is not None:
                if self.dtype != values.dtype:
                    raise ValueError(f"The provided dtype {dtype} doesn't match"
                                     f" the dtype of the provided values {values.dtype}")
            else:
                # Take dtype from the values
                self._dtype = values.dtype
            # Get the data
            if copy:
                values = values.copy()
            self._data = values._data

        elif len(values) == 0:
            # Return an empty Genotype Array
            if self._dtype is not None:
                self._data = np.array(values, dtype=GenotypeDtype._record_type)
            else:
                raise ValueError("Cannot create a Genotype Array with neither values nor a dtype")

        elif all([type(i) == Genotype for i in values]):
            # Sequence of Genotype objects
            genotype_array = self._from_sequence(scalars=values, dtype=dtype, copy=copy)
            # Replace self with the created array
            self._data = genotype_array._data
            self._dtype = genotype_array._dtype

        elif all([type(i) == str for i in values]):
            # List of Strings
            genotype_array = self._from_sequence_of_strings(strings=values, dtype=dtype, copy=copy)
            # Replace self with the created array
            self._data = genotype_array._data
            self._dtype = genotype_array._dtype

        else:
            raise ValueError(f"Unsupported `values` type passed to __init__: {type(values)}")

    # ExtensionArray Initialization Methods
    # -------------------------------------

    @classmethod
    def _from_sequence(cls, scalars, dtype: Optional[GenotypeDtype] = None, copy: bool = False) -> 'GenotypeArray':
        """
        Construct a new GenotypeArray from a sequence of Genotypes.
        Parameters
        ----------
        scalars : Sequence of Genotype objects
        dtype : dtype, optional
            Optional GenotypeDtype that must be compatible with the Genotypes
        copy : boolean, default False
            If True, copy the underlying data.
        
        Returns
        -------
        GenotypeArray
        """
        if len(scalars) == 0:
            # Pass the scalars object anyway, in case it is an empty GenotypeArray
            return cls(values=scalars, dtype=dtype)

        if dtype is None:
            # Use variant from first genotype
            variant = scalars[0].variant
        else:
            # Use the dtype variant
            variant = dtype.variant
        values = []
        for idx, gt in enumerate(scalars):
            if not variant.is_same_variant(gt.variant):
                raise ValueError(f"Variant for Genotype {idx} of {len(scalars)} ({gt.variant}) "
                                 f"is not compatible with the prior ones ({variant})")
            else:
                values.append((gt.allele1, gt.allele2))
        result = cls(values=[], dtype=GenotypeDtype(variant))
        result._data = np.array(values, dtype=GenotypeDtype._record_type)
        return result

    @classmethod
    def _from_sequence_of_strings(cls, strings, dtype, copy: bool = False):
        """Construct a new ExtensionArray from a sequence of strings.
        .. versionadded:: 0.24.0
        Parameters
        ----------
        strings : Sequence
            A list of genotype strings, like "A/C", "G/delG", etc
        dtype : dtype, optional
            GenotypeDtype with variant information used to process the strings
        copy : boolean, default False
            If True, copy the underlying data.
        Returns
        -------
        GenotypeArray
        """
        variant = dtype.variant
        return cls._from_sequence([variant.make_genotype_from_str(s) for s in strings], dtype, copy)

    @classmethod
    def _from_factorized(cls, values, original):
        """
        Reconstruct an ExtensionArray after factorization.
        Parameters
        ----------
        values : ndarray
            An integer ndarray with the factorized values.
        original : ExtensionArray
            The original ExtensionArray that factorize was called on.
        See Also
        --------
        factorize
        ExtensionArray.factorize
        """
        return cls(values, dtype=original.dtype)

    @classmethod
    def is_genotype_array(cls, other):
        if isinstance(other, cls):
            return True
        else:
            return False

    # Attributes
    # ----------
    @property
    def dtype(self) -> GenotypeDtype:
        """The specific parametized type"""
        return self._dtype

    @property
    def nbytes(self) -> int:
        """How many bytes to store this object in memory (2 per genotype)"""
        return 2 * len(self)

    def __getitem__(self, index):
        """
        Select a subset of self.
        Parameters
        ----------
        index : int, slice, or ndarray
            * int: The position in 'self' to get.
            * slice: A slice object, where 'start', 'stop', and 'step' are
              integers or None
            * ndarray: A 1-d boolean NumPy ndarray the same length as 'self'

        Returns
        -------
        item : Genotype or GenotypeArray

        Notes
        -----
        """
        # Check and convert the index
        index = pd.api.indexers.check_array_indexer(self._data, index)

        result = operator.getitem(self._data, index)

        if isinstance(result, np.ndarray):
            return GenotypeArray(values=result, dtype=self.dtype)
        elif isinstance(result, np.void):
            return Genotype(variant=self.dtype.variant, allele1=result['allele1'], allele2=result['allele2'])
        else:
            raise TypeError("Indexing error- unexpected type")

    def __setitem__(self, key: Union[int, np.ndarray], value: Union[Genotype, 'GenotypeArray']) -> None:
        if isinstance(value, Genotype):
            self._data[key] = (value.allele1, value.allele2)
        elif isinstance(value, GenotypeArray):
            self._data['allele1'][key] = value._data['allele1']
            self._data['allele2'][key] = value._data['allele2']
        elif isinstance(value, pd.Series) and isinstance(value.values, GenotypeArray):
            self._data['allele1'][key] = value.values._data['allele1']
            self._data['allele2'][key] = value.values._data['allele2']
        else:
            raise ValueError(f"Can't set the value in a GenotypeArray with '{type(value)}")

    def __len__(self):
        return len(self._data)

    def take(self, indices, allow_fill=False, fill_value=None):
        indices = np.asarray(indices, dtype='int')

        if allow_fill:
            # Unpack the genotype or use the default value for None
            fill_value = self.dtype.unpack_genotype(genotype=fill_value)

        if allow_fill:
            mask = (indices == -1)
            if len(self) == 0:
                if not (indices == -1).all():
                    raise IndexError("Invalid take for empty array. Must be all -1.")
                else:
                    print(indices)
                    took = (np.full((len(indices), 2), fill_value, dtype=self.dtype._record_type).reshape(-1))
                    return GenotypeArray(values=took, dtype=self.dtype)
            if (indices < -1).any():
                raise ValueError("Invalid value in 'indicies'. Must be all >= -1 for 'allow_fill=True'")

        took = self._data.take(indices)
        if allow_fill:
            took[mask] = fill_value

        return GenotypeArray(values=took, dtype=self.dtype)

    def copy(self):
        return GenotypeArray(self._data.copy(), self.dtype)

    def factorize(self, na_sentinel: int = -1) -> Tuple[np.ndarray, 'GenotypeArray']:
        """
        Return an array of ints indexing unique values
        """
        if len(self) == 0:
            return np.array([], dtype=np.int64), self

        codes = np.ones(len(self), dtype=np.int64)

        # Get list of unique genotypes in order they appear
        _, idx = np.unique(self._data, return_index=True)
        uniques = self._data[np.sort(idx)]
        uniques = GenotypeArray(values=uniques, dtype=self.dtype)

        # Update codes for unique values, not including NA
        for idx, gt in enumerate([gt for gt in uniques if gt != self.dtype.na_value]):
            codes[self == gt] = idx

        # Update codes for NA values
        codes[self.isna()] = na_sentinel

        # Return the codes and unique values (not including NA)
        return codes, uniques[~uniques.isna()]

    def unique(self) -> 'GenotypeArray':
        """Return a GenotypeArray of unique values"""
        unique_data = np.unique(self._data)
        return GenotypeArray(values=unique_data, dtype=self.dtype)

    def value_counts(self, dropna=True):
        raise NotImplementedError()

    #def _values_for_factorize(self):
    #    return self.astype(object), self.variant.make_genotype()

    def astype(self, dtype, copy=True):
        if isinstance(dtype, GenotypeDtype) or isinstance(dtype, object):
            if copy:
                return self.copy()
            else:
                return self
        else:
            raise ValueError(f"Can't coerce GenotypeArray to 'dtype'")

    def isna(self):
        """
        A 1-D array indicating if each value is missing
        """
        return (self._data['allele1'] == 255) & (self._data['allele2'] == 255)

    @classmethod
    def _concat_same_type(cls, to_concat):
        """
        Concatenate multiple array
        Parameters
        ----------
        to_concat : sequence of this type
        Returns
        -------
        ExtensionArray
        """
        # Check dtypes
        dtype = to_concat[0].dtype
        for a in to_concat:
            if a.dtype != dtype:
                raise ValueError(f"Incompatible types: {dtype} and {a.dtype}")

        data = np.concatenate([ga._data for ga in to_concat])

        return GenotypeArray(data, dtype)

    # Properties for the type parameters
    # ----------------------------------

    @property
    def variant(self):
        """
        Return the variant identifier
        """
        return self._dtype.variant

    # Operations
    # Note: genotypes are compared by first allele then second, using the order of alleles in the variant
    # ----------
    def _get_alleles_for_ops(self, other):
        """Get the scalar allele values (Genotype) or arrays of allele values (GenotypeArray)"""
        # Ensure the comparison is using the same variant
        if other.variant != self.variant:
            raise ValueError("Can't compare genotypes from different variants")
        # Get the alleles
        if isinstance(other, Genotype):
            # Get scalar values for alleles
            allele1 = other.allele1
            allele2 = other.allele2
        elif isinstance(other, GenotypeArray):
            # Get array values for alleles
            allele1 = other._data['allele1']
            allele2 = other._data['allele2']
        else:
            raise NotImplementedError

        return allele1, allele2

    def __eq__(self, other):
        allele1, allele2 = self._get_alleles_for_ops(other)
        return (self._data['allele1'] == allele1) & (self._data['allele2'] == allele2)

    def __lt__(self, other):
        allele1, allele2 = self._get_alleles_for_ops(other)
        return (self._data['allele1'] <= allele1) & (self._data['allele2'] < allele2)

    def __le__(self, other):
        allele1, allele2 = self._get_alleles_for_ops(other)
        return (self._data['allele1'] <= allele1) & (self._data['allele2'] <= allele2)

    def __gt__(self, other):
        allele1, allele2 = self._get_alleles_for_ops(other)
        return (self._data['allele1'] >= allele1) & (self._data['allele2'] > allele2)

    def __ge__(self, other):
        allele1, allele2 = self._get_alleles_for_ops(other)
        return (self._data['allele1'] >= allele1) & (self._data['allele2'] >= allele2)
