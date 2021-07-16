import operator
import re
from copy import copy
from typing import Dict, MutableMapping, Any, Optional, List, Union, Tuple, Iterable

import numpy as np
import pandas as pd
from pandas.core.arrays import ExtensionArray, BooleanArray, IntegerArray
from pandas.core.dtypes.dtypes import register_extension_dtype, PandasExtensionDtype
from pandas.core.dtypes.inference import is_list_like

from pandas_genomics.arrays.encoding_mixin import EncodingMixin
from pandas_genomics.arrays.info_mixin import InfoMixin
from pandas_genomics.scalars import Variant, Genotype, MISSING_IDX


@register_extension_dtype
class GenotypeDtype(PandasExtensionDtype):
    """
    An ExtensionDtype for genotype data.

    Parameters
    ----------
    variant: Variant or None
        The ~Variant associated with the genotype.
        If None, use an anonymous variant

    Attributes
    ----------
    variant: Variant
        The variant that the datatype is specific to

    Examples
    --------
    v = Variant(chromosome='12', position=112161652, id='rs12462', ref='T', alt=['C',], score=25)
    >>> GenotypeDtype(v)
    genotype(2n)[12; 112161652; rs12462; T; C]Q25
    """

    # Internal attributes
    # -------------------
    # metadata field names
    _metadata = ("variant",)
    # Regular expression
    # TODO: Validate ref/alt more specifically
    _match = re.compile(
        r"genotype\((?P<ploidy>[0-9]+)n\)\["
        r"(?P<chromosome>.+); "
        r"(?P<position>[0-9]+); "
        r"(?P<id>.+); "
        r"(?P<ref>.+); "
        r"(?P<alt>.+)\]"
        r"(Q(?P<score>[0-9]+))?$"
    )
    kind = "O"
    type = Genotype

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
    def __init__(self, variant: Optional[Variant] = None):
        # Set variant
        if variant is None:
            variant = Variant()
        self.variant = variant

        # Data backing the GenotypeArray is stored as a numpy structured array
        # An unsigned integer for each allele in the genotype indexing the list of possible alleles
        # An unsigned integer for the genotype score (255 if missing)
        self._record_type = np.dtype(
            [
                ("allele_idxs", np.uint8, (self.variant.ploidy,)),
                ("gt_score", np.uint8),
            ]
        )
        self.itemsize = self._record_type.itemsize

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
            Should be formatted like `genotype(<ploidy>n)[<chromosome>; <position>; <id>; <ref>; <alt>]`

        Examples
        --------
        >>> GenotypeDtype.construct_from_string('genotype(2n)[chr1; 123456; rs12345; A; T,G]')
        genotype(2n)[chr1; 123456; rs12345; A; T,G]
        """
        if isinstance(string, str):
            msg = "Cannot construct a 'GenotypeDtype' from '{}'"
            try:
                match = cls._match.match(string)
                if match is not None:
                    d = match.groupdict()
                    # Score is optional, so it may be None
                    score = d["score"]
                    if score is not None:
                        score = int(score)
                    variant = Variant(
                        chromosome=d["chromosome"],
                        position=int(d["position"]),
                        id=d["id"],
                        ref=d["ref"],
                        alt=d["alt"].split(","),
                        ploidy=int(d["ploidy"]),
                        score=score,
                    )
                    return cls(variant=variant)
                else:
                    raise TypeError(msg.format(string))
            except Exception:
                raise TypeError(msg.format(string))
        else:
            raise TypeError(
                f"'construct_from_string' expects a string, got {type(string)}>"
            )

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
        >>> genotype = variant.make_genotype_from_str('C/T', add_alleles=True)
        >>> GenotypeDtype.from_genotype(genotype)
        genotype(2n)[12; 112161652; rs12462; ref=N; alt=T,C]
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
            if dtype.lower().startswith("genotype("):
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
        variant_score_str = ""
        if self.variant.score is not None:
            variant_score_str = f"Q{self.variant.score}"
        return (
            f"genotype({self.variant.ploidy}n)["
            f"{self.variant.chromosome}; "
            f"{self.variant.position}; "
            f"{self.variant.id}; "
            f"{self.variant.ref}; "
            f"{self.variant.alt}]" + variant_score_str
        )

    def __hash__(self):
        return hash(str(self))

    def __copy__(self):
        """Create a copy to avoid references to the same variant"""
        return GenotypeDtype(copy(self.variant))

    def __eq__(self, other):
        if isinstance(other, str):
            return other == self.name

        return isinstance(other, GenotypeDtype) and self.variant == other.variant

    # Pickle compatibility
    # --------------------
    def __getstate__(self) -> Dict[str, Any]:
        # pickle support; we don't want to pickle the cache
        return {k: getattr(self, k, None) for k in self._metadata}

    def __setstate__(self, state: MutableMapping[str, Any]) -> None:
        self.variant = state.pop("variant")

    # Other internal methods
    # ----------------------

    def unpack_genotype(self, genotype: Optional[Genotype]):
        if genotype is None:
            return self.unpack_genotype(self.na_value)
        elif not isinstance(genotype, Genotype):
            raise ValueError(
                f"The fill value must be None or a Genotype object, not {type(genotype)}"
            )
        else:
            score = float("NaN")
            if genotype.score is not None:
                score = genotype.score
            return np.array(
                list(genotype.allele_idxs) + [score], dtype=self._record_type
            )


class GenotypeArray(ExtensionArray, EncodingMixin, InfoMixin):
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
    data: np.dtype("u8") with shape (<genotypes>, <ploidy>)
        The genotype values encoded as indices into the allele list of the dtype
    """

    # array priority higher than numpy scalars
    __array_priority__ = 1000

    def __init__(
        self,
        values: Union[List[Genotype], "GenotypeArray", np.ndarray],
        dtype: Optional[GenotypeDtype] = None,
        copy: bool = False,
    ):
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
        if not is_list_like(values):
            values = [values]

        if isinstance(values, np.ndarray) and (
            values.dtype == self._dtype._record_type
        ):
            # Stored data format
            self._data = values

        elif self.is_genotype_array(values):
            # values is a GenotypeArray, simply check the dtype and return
            if self.dtype is not None:
                if self.dtype != values.dtype:
                    raise ValueError(
                        f"The provided dtype {dtype} doesn't match"
                        f" the dtype of the provided values {values.dtype}"
                    )
            else:
                # Take dtype from the values
                self._dtype = values.dtype
            # Get the data
            if copy:
                values = values.copy()
            self._data = values._data

        elif len(values) == 0:
            # Return an empty Genotype Array with the given GenotypeDtype
            self._data = np.array(values, dtype=self._dtype._record_type)

        elif all([type(i) == Genotype for i in values]):
            # Sequence of Genotype objects
            genotype_array = self._from_sequence(scalars=values, dtype=dtype, copy=copy)
            # Replace self with the created array
            self._data = genotype_array._data
            self._dtype = genotype_array._dtype

        elif all([type(i) == str for i in values]):
            # List of Strings
            genotype_array = self._from_sequence_of_strings(
                strings=values, dtype=dtype, copy=copy
            )
            # Replace self with the created array
            self._data = genotype_array._data
            self._dtype = genotype_array._dtype

        else:
            raise ValueError(
                f"Unsupported `values` type passed to __init__: {type(values)}"
            )

        # Set an anonymous dtype if one was not set
        if self.dtype is None:
            self._dtype = GenotypeDtype()

    # ExtensionArray Initialization Methods
    # -------------------------------------

    @classmethod
    def _from_sequence(
        cls,
        scalars: Union[Genotype, Iterable[Genotype]],
        dtype: Optional[GenotypeDtype] = None,
        copy: bool = False,
    ) -> "GenotypeArray":
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
        if type(scalars) == Genotype:
            # Despite the documentation, scalars is not always a sequence of objects, sometimes it is just one
            scalars = [scalars]

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
            if not variant.is_same_position(gt.variant):
                raise ValueError(
                    f"Variant for Genotype {idx} of {len(scalars)} ({gt.variant}) "
                    f"is not compatible with the prior ones ({variant})"
                )
            elif variant.score != gt.variant.score:
                raise ValueError(
                    f"Variant for Genotype {idx} of {len(scalars)} ({gt.variant}) "
                    f"is compatible, but has a different variant score"
                )
            else:
                values.append((gt.allele_idxs, gt._float_score))
        result = cls(values=[], dtype=GenotypeDtype(variant))
        result._data = np.array(values, dtype=result._dtype._record_type)
        return result

    @classmethod
    def _from_sequence_of_strings(cls, strings, dtype, copy: bool = False):
        """Construct a new ExtensionArray from a sequence of strings.

        Parameters
        ----------
        strings : Sequence
            A list of genotype strings, like "A/C", "G/delG", etc
        dtype : dtype, optional
            GenotypeDtype with variant information used to process the strings
        copy : boolean, default False
            If True, copy the underlying data.

        Note
        ----
        This will not allow the use of genotype scores as they are not encoded in genotype strings

        Returns
        -------
        GenotypeArray
        """
        variant = dtype.variant
        return cls._from_sequence(
            [variant.make_genotype_from_str(s) for s in strings], dtype, copy
        )

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
        """How many bytes to store this object in memory"""
        return self._dtype.itemsize * len(self)

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
            variant = self.dtype.variant
            return Genotype(
                variant=variant,
                allele_idxs=result["allele_idxs"],
                score=result["gt_score"] if not np.isnan(result["gt_score"]) else None,
            )
        else:
            raise TypeError("Indexing error- unexpected type")

    def __setitem__(
        self,
        key: Union[int, np.ndarray],
        value: Union[Genotype, "GenotypeArray", List[Genotype]],
    ) -> None:
        if isinstance(value, list):
            # Convert list to genotype array, throwing an error if it doesn't work
            value = self._from_sequence(value, dtype=self.dtype)

        # Validate the key
        if isinstance(key, List):
            key = pd.Series(key)
            if key.isna().sum() > 0:
                raise ValueError(
                    "Cannot index with an integer indexer containing NA values"
                )
        if isinstance(key, BooleanArray):
            # Convert to a normal boolean array after making NaN rows False
            key = key.fillna(False).astype("bool")
        # Handle pandas IntegerArray
        if isinstance(key, IntegerArray):
            if key.isna().sum() > 0:
                # Raise an error if there are any NA values
                raise ValueError(
                    "Cannot index with an integer indexer containing NA values"
                )
            else:
                # Convert to a regular numpy array of ints
                key = key.astype("int")
        # Ensure a mask doesn't have an incorrect length
        if isinstance(key, np.ndarray) and key.dtype == "bool":
            if len(key) != len(self):
                raise IndexError("wrong length")
        # Update allele values directly
        if isinstance(value, Genotype):
            self._data[key] = tuple([value.allele_idxs, value._float_score])
        elif isinstance(value, GenotypeArray):
            self._data[key] = value._data
        elif isinstance(value, pd.Series) and isinstance(value.values, GenotypeArray):
            self._data[key] = value.values._data
        else:
            raise ValueError(
                f"Can't set the value in a GenotypeArray with '{type(value)}"
            )

    def __len__(self):
        return len(self._data)

    def take(self, indexer, allow_fill=False, fill_value=None):
        indexer = np.asarray(indexer)
        msg = (
            "Index is out of bounds or cannot do a "
            "non-empty take from an empty array."
        )

        if allow_fill:
            if fill_value is None:
                fill_value = self.dtype.na_value
            # bounds check
            if (indexer < -1).any():
                raise ValueError
            try:
                output = [self[loc] if loc != -1 else fill_value for loc in indexer]
            except IndexError as err:
                raise IndexError(msg) from err
        else:
            try:
                output = [self[loc] for loc in indexer]
            except IndexError as err:
                raise IndexError(msg) from err

        return self._from_sequence(scalars=output, dtype=self.dtype)

    def copy(self):
        return GenotypeArray(self._data.copy(), copy(self.dtype))

    def factorize(self, na_sentinel: int = -1) -> Tuple[np.ndarray, "GenotypeArray"]:
        """
        Return an array of ints indexing unique values
        """
        if len(self) == 0:
            return np.array([], dtype=np.int64), self

        codes = np.ones(len(self), dtype=np.int64)

        # Get list of unique genotypes in order they appear (not counting scores)
        _, idx = np.unique(self.allele_idxs, return_index=True, axis=0)
        uniques = self._data[np.sort(idx)]
        uniques = GenotypeArray(values=uniques, dtype=self.dtype)

        # Update codes for unique values, not including NA
        for idx, gt in enumerate([gt for gt in uniques if gt != self.dtype.na_value]):
            codes[self == gt] = idx

        # Update codes for NA values
        codes[self.isna()] = na_sentinel

        # Return the codes and unique values (not including NA)
        return codes, uniques[~uniques.isna()]

    def unique(self) -> "GenotypeArray":
        """Return a GenotypeArray of unique values"""
        _, idx = np.unique(self.allele_idxs, return_index=True, axis=0)
        return GenotypeArray(values=self._data[idx], dtype=self.dtype)

    def value_counts(self, dropna=True):
        """Return a Series of unique counts with a GenotypeArray index"""
        _, unique_idx, unique_counts = np.unique(
            self.allele_idxs, return_index=True, return_counts=True, axis=0
        )
        result = pd.Series(
            unique_counts,
            index=GenotypeArray(values=self._data[unique_idx], dtype=self.dtype),
        )
        if dropna:
            result = result.loc[result.index != self.dtype.na_value]
        return result

    def astype(self, dtype, copy=True):
        if isinstance(dtype, GenotypeDtype):
            if copy:
                self = self.copy()
            return self
        return super(GenotypeArray, self).astype(dtype)

    def isna(self):
        """
        A 1-D array indicating if each value is missing
        """
        return (self.allele_idxs == MISSING_IDX).all(axis=1)

    @classmethod
    def _concat_same_type(cls, to_concat, axis: int = 0):
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
        dtypes = {a.dtype for a in to_concat}
        if len(dtypes) != 1:
            raise ValueError(
                "to_concat must have the same dtype for all values", dtypes
            )

        data = np.concatenate([ga._data for ga in to_concat], axis=axis)

        return GenotypeArray(data, list(dtypes)[0])

    # Properties for accessing array metadata
    # ----------------------------------------
    @property
    def variant(self):
        """
        Return the variant identifier
        """
        return self._dtype.variant

    # Properties for accessing individual portions of the internal _data array
    # ------------------------------------------------------------------------
    @property
    def allele_idxs(self):
        """
        Return the allele indices for each genotype
        """
        return self._data["allele_idxs"]

    @property
    def gt_scores(self):
        """
        Return the genotype score for each genotype (as a float)
        """
        scores = self._data["gt_score"].copy().astype("float")
        scores[scores == MISSING_IDX] = np.nan
        return scores

    # Operations
    # Note: genotypes are compared by first allele then second, using the order of alleles in the variant
    # ----------
    def _get_alleles_for_ops(self, other):
        """
        Get the scalar allele values (Genotype or str)
                or arrays of allele values (GenotypeArray)
                or None (NotImplemented)
        """
        if isinstance(other, Genotype) or isinstance(other, GenotypeArray):
            # Get scalar values for alleles
            allele_idxs = other.allele_idxs
            # Ensure the comparison is using the same variant
            if self.variant != other.variant:
                return None
        elif type(other) == str:
            allele_idxs = tuple(
                self.variant.get_idx_from_allele(a) for a in other.split("/")
            )
        else:
            return None
        return allele_idxs

    def __eq__(self, other):
        allele_idxs = self._get_alleles_for_ops(other)
        if allele_idxs is None:
            return NotImplemented
        return (self.allele_idxs == allele_idxs).all(axis=1)

    def __ne__(self, other):
        allele_idxs = self._get_alleles_for_ops(other)
        if allele_idxs is None:
            return NotImplemented
        return (self.allele_idxs != allele_idxs).any(axis=1)

    def __lt__(self, other):
        allele_idxs = self._get_alleles_for_ops(other)
        if allele_idxs is None:
            return NotImplemented
        return (self.allele_idxs <= allele_idxs).all(axis=1) & (
            self.allele_idxs < allele_idxs
        ).any(axis=1)

    def __le__(self, other):
        allele_idxs = self._get_alleles_for_ops(other)
        if allele_idxs is None:
            return NotImplemented
        return (self.allele_idxs <= allele_idxs).all(axis=1)

    def __gt__(self, other):
        allele_idxs = self._get_alleles_for_ops(other)
        if allele_idxs is None:
            return NotImplemented
        return (self.allele_idxs >= allele_idxs).all(axis=1) & (
            self.allele_idxs > allele_idxs
        ).any(axis=1)

    def __ge__(self, other):
        allele_idxs = self._get_alleles_for_ops(other)
        if allele_idxs is None:
            return NotImplemented
        return (self.allele_idxs >= allele_idxs).all(axis=1)

    ###################
    # Utility Methods #
    ###################
    def set_reference(self, allele: Union[str, int]) -> None:
        """
        Change the reference allele (in-place) by specifying an allele index value or an allele string

        Parameters
        ----------
        allele: int or str
            The allele that will be set as the reference allele.
            Either the allele string, or the index into the variant allele list

        Returns
        -------
        None
        """
        # Get the allele as an integer and as a string
        if type(allele) == str:
            allele_idx = self.variant.get_idx_from_allele(allele, add=False)
            allele_str = allele
        elif type(allele) == int:
            if not self.variant.is_valid_allele_idx(allele):
                raise ValueError(
                    f"{allele} is not a valid allele index,"
                    f" the variant has {len(self.variant.alleles)} alleles."
                )
            allele_idx = allele
            allele_str = self.variant.alleles[allele]
        else:
            raise ValueError(
                f"The `allele` must be a str or int, not an instance of '{type(allele)}'"
            )

        if allele_idx == 0:
            # Nothing to do, this is already the reference
            return self

        # Update the list of alleles
        old_ref = self.variant.alleles[0]
        # Replace existing value with the old ref
        self.variant.alleles[allele_idx] = old_ref
        # Add new ref to the beginning and remove the old ref
        self.variant.alleles = [
            allele_str,
        ] + self.variant.alleles[1:]

        # Update stored alleles
        was_ref = self._data["allele_idxs"] == 0
        was_allele = self._data["allele_idxs"] == allele_idx
        # What was the reference is now the new reference position
        self._data["allele_idxs"][was_ref] = allele_idx
        # What was the allele is now reference (0)
        self._data["allele_idxs"][was_allele] = 0
