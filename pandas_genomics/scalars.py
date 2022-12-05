"""
Scalars
-------
This module contains scalar types, some of which are used in the ExtensionArrays.  They may also be useful on their own.

.. autosummary::
     :toctree: scalars

     Variant
     Genotype
     Region
"""
import uuid
from dataclasses import dataclass
from typing import List, Optional, Tuple, Union

# Integer indicating a missing allele or genotype score.
# Each variant must have 254 alleles max and the maximum genotype score is 254.
MISSING_IDX = 255


class Variant:
    """
    Information about a variant.

    Parameters
    ----------
    chromosome: str, optional
        None by default, through this usually is not desired.
    position: int, optional
        (1-based, the default is 0 for none/unknown)
    id: str, optional
        None by default
    ref: str, optional
        Reference allele, 'N' by default
    alt: List[str] or str, optional
        One (or a list of) possible alternate alleles, empty list by default
    ploidy: int, optional
        The number of alleles that should be present.  Assumed to be 2 if not specified.
    score: int, optional
        A quality score for the Variant.  No assumptions are made about the meaning.

    Notes
    -----
    Missing alleles are shown as `.`, such as `A/.` or `./.`


    Examples
    --------
    >>> variant = Variant('12', 112161652, id='rs12462', ref='A', alt=['C', 'T'])
    >>> print(variant)
    rs12462[chr=12;pos=112161652;ref=A;alt=C,T]
    """

    def __init__(
        self,
        chromosome: Optional[str] = None,
        position: int = 0,
        id: Optional[str] = None,
        ref: Optional[str] = "N",
        alt: Optional[Union[List[str], str]] = None,
        ploidy: Optional[int] = None,
        score: Optional[int] = None,
    ):
        self.chromosome = chromosome
        self.position = position
        if id is None:
            # Use a UUID to avoid duplicate IDs
            id = str(uuid.uuid4())
        self.id = id
        if ploidy is None:
            self.ploidy = 2
        elif ploidy < 1:
            raise ValueError(
                f"If specified, ploidy must be >= 1, but was specified as {ploidy}"
            )
        else:
            self.ploidy = ploidy
        self.score = score

        # Validate alleles
        # TODO: Validate ref and alt alleles using regex
        if ref is None:
            ref = "N"
        if alt is None:
            alt = []
        elif type(alt) == str:
            alt = [
                alt,
            ]
        if ref in alt:
            raise ValueError(
                f"The ref allele ({ref}) was also listed as an alt allele."
            )
        if len(alt) >= MISSING_IDX:
            raise ValueError(
                f"Too many alternate alleles ({len(alt):,} > {MISSING_IDX})"
            )

        # Store alleles in a big list with the ref first
        self.alleles = [
            ref,
        ] + alt

        # Validate the passed parameters
        if self.chromosome is not None and (
            ";" in self.chromosome or "," in self.chromosome
        ):
            raise ValueError(
                f"The chromosome cannot contain ';' or ',': '{self.chromosome}'"
            )
        if self.position > ((2**31) - 2):
            raise ValueError(
                f"The 'position' value may not exceed 2^31-2, {self.position:,} was specified"
            )
        if self.id is not None and (";" in self.id or "," in self.id):
            raise ValueError(f"The variant_id cannot contain ';' or ',': '{self.id}'")

    @property
    def ref(self):
        return self.alleles[0]

    @property
    def alt(self):
        return ",".join(self.alleles[1:])

    def __str__(self):
        score_str = ""
        if self.score is not None:
            score_str = f"Q{self.score}"
        return f"{self.id}[chr={self.chromosome};pos={self.position};ref={self.ref};alt={self.alt}]{score_str}"

    def __repr__(self):
        return (
            f"Variant(chromosome={self.chromosome}, position={self.position}, "
            f"id={self.id}, ref={self.ref}, alt={self.alleles[1:]}, "
            f"ploidy={self.ploidy}, score={self.score})"
        )

    def __eq__(self, other):
        if other.__class__ is not self.__class__:
            return NotImplemented
        return (
            (self.chromosome == other.chromosome)
            & (self.position == other.position)
            & (self.id == other.id)
            & (self.alleles == other.alleles)
            & (self.ploidy == other.ploidy)
            & (self.score == other.score)
        )

    def add_allele(self, allele):
        """
        Add a potential allele to the variant

        Parameters
        ----------
        allele: str

        Returns
        -------
        None

        """
        if allele in self.alleles:
            raise ValueError("Allele already exists in the variant")
        if len(self.alleles) < MISSING_IDX:
            self.alleles.append(allele)
        else:
            raise ValueError(
                f"Couldn't add new allele to {self}, {MISSING_IDX} alleles max."
            )
        print(f"Alternate alleles = {self.alt}")

    def get_idx_from_allele(self, allele: Optional[str], add: bool = False) -> int:
        """
        Get the integer value for an allele based on this variant's list of alleles,
        optionally adding it as a new allele

        Parameters
        ----------
        allele: str
        add: bool
            If `add` is True and the allele doesn't exist in the variant, add it
            If 'add' is False and the allele doesn't exist in the variant, return the missing index value

        Returns
        -------
        int
            The integer value for the allele in the variant

        """
        if allele is None or (allele == "."):
            return MISSING_IDX
        else:
            try:
                # Get allele idx
                allele_idx = self.alleles.index(allele)
            except ValueError:
                if add:
                    # Add as a new allele
                    self.add_allele(allele)
                    allele_idx = len(self.alleles) - 1
                else:
                    raise ValueError(f"'{allele}' is not an allele in {self}")
            return allele_idx

    def get_allele_from_idx(self, idx: int) -> str:
        """
        Get the allele corresponding to an index value for this variant

        Parameters
        ----------
        idx: int

        Returns
        -------
        str
            String representation of the allele
        """

        if idx == MISSING_IDX:
            return "."  # Missing allele symbol, same as VCF
        elif idx > len(self.alleles) - 1:
            raise ValueError(
                f"Invalid index (idx) for a variant with {len(self.alleles)} alleles"
            )
        else:
            return self.alleles[idx]

    def is_valid_allele_idx(self, idx: int) -> bool:
        """
        Validate the integer value for an allele with respect to this variant

        Parameters
        ----------
        idx: int
            allele index value to be checked

        Returns
        -------
        bool
            True if the index value is valid (including None/Missing value)
            False if the index value is not valid (negative or a nonexistent allele index)

        """
        if idx == MISSING_IDX:
            # None/Missing
            return True
        elif idx < 0:
            return False
        elif idx > (len(self.alleles) - 1):
            return False
        else:
            return True

    def is_same_position(self, other):
        """
        Confirms this is a variant at the same position

        Parameters
        ----------
        other: Variant
            Another variant to compare

        Returns
        -------
        bool
            True if the variant_id, chromosome, position, and ref allele all match.
            False otherwise.

        """
        if isinstance(other, Variant):
            return (
                self.id == other.id
                and self.chromosome == other.chromosome
                and self.position == other.position
                and self.alleles[0] == other.alleles[0]
                and self.ploidy == other.ploidy
            )
        else:
            return False

    def make_genotype(self, *alleles: str, add_alleles: bool = False) -> "Genotype":
        """
        Create a Genotype object associated with this variant using the specified allele(s)

        Parameters
        ----------
        alleles:
            zero or more alleles (as strings) that make up the genotype.
            If ploidy > the number of specified alleles, missing alleles will be filled in.
        add_alleles: bool
            False by default.  If True, add alleles if the Variant doesn't have them yet.

        Returns
        -------
        Genotype
            A Genotype based on this variant with the specified alleles
        """
        if len(alleles) > self.ploidy:
            raise ValueError(
                f"Too many alleles ({len(alleles)} specified for a variant with ploidy of {self.ploidy}"
            )
        allele_idxs = [self.get_idx_from_allele(a, add=add_alleles) for a in alleles]
        missing_idxs = [MISSING_IDX] * (self.ploidy - len(self.alleles))
        return Genotype(self, allele_idxs + missing_idxs)

    def make_genotype_from_str(
        self, gt_str: str, sep: str = "/", add_alleles: bool = False
    ) -> "Genotype":
        """
        Create a Genotype object associated with this variant using the specified allele string

        Parameters
        ----------
        gt_str: str
            A string containing two alleles separated by a delimiter character
        sep: str
            The delimiter used to separate alleles in the gt_str ('/' by default)
        add_alleles: bool
            False by default.  If True, add alleles if the Variant doesn't have them yet.

        Returns
        -------
        Genotype
            A Genotype based on this variant with the specified alleles
        """
        # Process Allele String
        alleles = gt_str.split(sep)
        if len(alleles) > self.ploidy:
            raise ValueError(
                f"Too many alleles ({len(alleles)} specified for a variant with ploidy of {self.ploidy}"
            )
        allele_idxs = [self.get_idx_from_allele(a, add=add_alleles) for a in alleles]
        missing_idxs = [MISSING_IDX] * (self.ploidy - len(self.alleles))
        return Genotype(self, allele_idxs + missing_idxs)

    def as_dict(self):
        """Return the variant information as a dictionary"""
        return {
            "id": self.id,
            "chromsome": self.chromosome,
            "position": self.position,
            "ref": self.ref,
            "alt": self.alt,
            "ploidy": self.ploidy,
            "score": self.score,
        }


class Genotype:
    """
    Genotype information associated with a specific variant.
    Defaults to using an anonymous variant with two unknown alleles (diploid).
    Usually created with methods on ~Variant

    Parameters
    ----------
    variant: pandas_genomics.scalars.variant.Variant
    allele_idxs: List[int]
        Alleles encoded as indexes into the variant allele list
    score: int, optional
        A quality score for the Genotype between 0 and 254.  255 or < 0 is treated as missing.

    Examples
    --------
    >>> variant = Variant('12', 112161652, 'rs12462')
    >>> genotype = variant.make_genotype_from_str('C/T')
    >>> print(genotype)
    C/T

    >>> missing_genotype = Genotype(variant)
    >>> print(missing_genotype)
    <Missing>
    """

    def __init__(
        self,
        variant: Variant,
        allele_idxs: Optional[Union[Tuple[int], List[int]]] = None,
        score: Optional[int] = None,
    ):

        # Determine alleles/ploidy
        if allele_idxs is not None and len(allele_idxs) > variant.ploidy:
            raise ValueError(
                f"{len(allele_idxs)} alleles given for a variant with ploidy of {variant.ploidy}"
            )
        elif allele_idxs is None:
            allele_idxs = []

        # Fill in any missing allele_idxs
        allele_idxs = list(allele_idxs) + [MISSING_IDX] * (
            variant.ploidy - len(allele_idxs)
        )

        # Ensure allele_idxs is a sorted tuple
        allele_idxs = tuple(sorted(allele_idxs))

        self.variant = variant
        self.allele_idxs = allele_idxs
        if score is not None:
            score = int(score)
            if score < 0 or score > 255:
                raise ValueError("The score must be between 0 and 255, inclusive")
            elif score == 255:
                score = None
        self.score = score

        # Validate parameters
        for a in self.allele_idxs:
            if not self.variant.is_valid_allele_idx(a):
                raise ValueError(f"Invalid allele index for {self.variant}: {a}")

    def __str__(self):
        if all([i == MISSING_IDX for i in self.allele_idxs]):
            return "<Missing>"
        else:
            return "/".join(
                [self.variant.get_allele_from_idx(i) for i in self.allele_idxs]
            )

    def __repr__(self):
        score_str = ""
        if self.score is not None:
            score_str = f"Q{self.score}"
        return f"Genotype(variant={self.variant})[{str(self)}]" + score_str

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, other):
        if other.__class__ is not self.__class__:
            return NotImplemented
        if self.variant != other.variant:
            raise NotImplementedError("Can't compare different variants")
        return self.allele_idxs == other.allele_idxs

    def __lt__(self, other):
        if other.__class__ is not self.__class__:
            return NotImplemented
        if self.variant != other.variant:
            raise NotImplementedError("Can't compare different variants")
        else:
            # Compare allele index values for sorting
            for selfa, othera in zip(self.allele_idxs, other.allele_idxs):
                if selfa < othera:
                    return True
                elif othera < selfa:
                    return False
            # All Equal
            return False

    def __gt__(self, other):
        if other.__class__ is not self.__class__:
            return NotImplemented
        if self.variant != other.variant:
            raise NotImplementedError("Can't compare different variants")
        else:
            # Compare allele index values for sorting
            for selfa, othera in zip(self.allele_idxs, other.allele_idxs):
                if selfa > othera:
                    return True
                elif othera > selfa:
                    return False
            # All Equal
            return False

    def __le__(self, other):
        if other.__class__ is not self.__class__:
            return NotImplemented
        if self.variant != other.variant:
            raise NotImplementedError("Can't compare different variants")
        else:
            # Compare allele index values for sorting
            for selfa, othera in zip(self.allele_idxs, other.allele_idxs):
                if selfa < othera:
                    return True
                elif othera < selfa:
                    return False
            # All Equal
            return True

    def __ge__(self, other):
        if other.__class__ is not self.__class__:
            return NotImplemented
        if self.variant != other.variant:
            raise NotImplementedError("Can't compare different variants")
        else:
            # Compare allele index values for sorting
            for selfa, othera in zip(self.allele_idxs, other.allele_idxs):
                if selfa > othera:
                    return True
                elif othera > selfa:
                    return False
            # All Equal
            return True

    def is_missing(self) -> bool:
        """
        Returns
        -------
        bool
            True if the variant is missing (both alleles are None), otherwise False
        """
        return all([a == MISSING_IDX for a in self.allele_idxs])

    @property
    def _float_score(self):
        """Convenience method for storing score as a uint8"""
        if self.score is None:
            return 255
        else:
            return self.score


@dataclass(order=True)
class Region:
    chromosome: str
    start: int
    end: int
    name: str = ""
    """
    Information associated with a genomic region.

    Parameters
    ----------
    chomosome: str
    start, end: int
        1-based, open-ended (includes start, not end)
    name: str
        name/comment for the region

    Examples
    --------
    >>> region1 = Region('chr12', 12345, 12387)
    >>> region2 = Region('chr12', 5236373, 5246380)
    >>> print(region1 < region2)
    True
    """

    def __post_init__(self):
        """Run after init, use to validate input"""
        if not isinstance(self.chromosome, str):
            raise TypeError(f"chromosome should be a str, not {type(self.chromosome)}")
        if not isinstance(self.start, int) or not isinstance(self.end, int):
            raise TypeError("start and end must be integers")
        if self.start < 1 or self.end < 1:
            raise ValueError("start and end must both be > 0 (1-based indexing)")
        if self.start >= self.end:
            raise ValueError(
                "start must be <= end, since the start is included and the end is not."
            )

    def contains_variant(self, var: Variant):
        if var.chromosome != self.chromosome:
            return False
        if var.position < self.start:
            return False
        if var.position >= self.end:
            return False
        return True
