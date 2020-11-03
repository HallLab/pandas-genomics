"""
Scalars
-------
This module contains scalar types used in the ExtensionArrays.  They may also be useful on their own.

.. autosummary::
     :toctree: scalars

     Variant
     Genotype
"""

from typing import Optional, List, Tuple, Union

MISSING_IDX = (
    255  # Integer indicating a missing allele.  Each variant must have 254 alleles max.
)


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
    alt: List[str], optional
        List of possible alternate alleles, empty list by default

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
        alt: Optional[List[str]] = None,
    ):
        self.chromosome = chromosome
        self.position = position
        self.id = id

        # Validate alleles
        # TODO: Validate ref and alt alleles using regex
        if ref is None:
            ref = "N"
        if alt is None:
            alt = []
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
        if self.position > ((2 ** 31) - 2):
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
        return f"{self.id}[chr={self.chromosome};pos={self.position};ref={self.ref};alt={self.alt}]"

    def __repr__(self):
        return (
            f"Variant(chromosome={self.chromosome}, position={self.position},"
            f"id={self.id}, ref={self.ref}, alt={self.alleles[1:]})"
        )

    def __eq__(self, other):
        if other.__class__ is not self.__class__:
            return NotImplemented
        return (
            (self.chromosome == other.chromosome)
            & (self.position == other.position)
            & (self.id == other.id)
            & (self.alleles == other.alleles)
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
            raise ValueError(f"Allele already exists in the variant")
        if len(self.alleles) < MISSING_IDX:
            self.alleles.append(allele)
        else:
            raise ValueError(
                f"Couldn't add new allele to {self}, {MISSING_IDX} alleles max."
            )
        print(f"Alternate alleles = {self.alt}")

    def get_allele_idx(self, allele: Optional[str], add: bool = False) -> int:
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
        if allele is None:
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

    def is_same_variant(self, other):
        """
        Confirms this is the same variant, other than the allele list.

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
            )
        else:
            return False

    def make_genotype(self, *alleles: str, add_alleles: bool = False) -> "Genotype":
        """
        Create a Genotype object associated with this variant using the specified allele(s)

        Parameters
        ----------
        alleles:
            one or more alleles (as strings) that make up the genotype
        add_alleles: bool
            False by default.  If True, add alleles if the Variant doesn't have them yet.

        Returns
        -------
        Genotype
            A Genotype based on this variant with the specified alleles
        """
        if len(alleles) == 0:
            # Create diploid missing genotype
            return Genotype(self, [MISSING_IDX, MISSING_IDX])
        else:
            allele_idxs = [self.get_allele_idx(a, add=add_alleles) for a in alleles]
            return Genotype(self, allele_idxs)

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
        allele_idxs = [self.get_allele_idx(a, add=add_alleles) for a in alleles]
        return Genotype(self, allele_idxs)

    def make_genotype_from_plink_bits(self, plink_bits: str) -> "Genotype":
        """
        Create a genotype from PLINK Bed file bits.
        Assumes the Variant has only the first and second alleles from the matching bim file

        Parameters
        ----------
        plink_bits: str
            A string with allele indices as encoded in plink format, one of {'00', '01', '10', '11'}

        Returns
        -------
        Genotype
            A Genotype based on this variant with the specified alleles
        """
        # Raise an error if the variant has more than a ref and alt allele
        if len(self.alleles) != 2:
            raise ValueError(
                "Genotypes can only be created from plink bitcodes if there are exactly two alleles"
            )
        # Process Allele String
        if plink_bits == "00":
            a1 = 0
            a2 = 0
        elif plink_bits == "01":
            a1 = MISSING_IDX
            a2 = MISSING_IDX
        elif plink_bits == "10":
            a1 = 0
            a2 = 1
        elif plink_bits == "11":
            a1 = 1
            a2 = 1
        else:
            raise ValueError(f"Invalid plink_bits: '{plink_bits}'")

        return Genotype(self, [a1, a2])

    def make_genotype_from_vcf_record(self, vcf_record: Tuple) -> "Genotype":
        """
        Create a genotype from VCF records loaded as cyvcf2.VCF().genotypes

        Parameters
        ----------
        vcf_record: Tuple
            The list is an array of allele indicies (where "-1" is missing) with a boolean for phased status at the nd

        Returns
        -------
        Genotype
            A Genotype based on this variant with the specified alleles
        """
        # TODO: Replace with pysam?
        allele_idxs = vcf_record[:-1]
        # is_phased = vcf_record[-1]

        # Replace -1 with MISSING_IDX
        allele_idxs = [a if a != -1 else MISSING_IDX for a in allele_idxs]

        for a in allele_idxs:
            assert self.is_valid_allele_idx(a)

        return Genotype(self, allele_idxs)


class Genotype:
    """
    Genotype information associated with a specific variant.
    Defaults to using an anonymous variant with two unknown alleles (diploid).
    Usually created with methods on ~Variant

    Parameters
    ----------
    variant: pandas_genomics.scalars.variant.Variant
    allele_idxs: Alleles encoded as indexes into the variant allele list
    ploidy: Number of alleles in the genotype

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
        ploidy: Optional[int] = None,
    ):

        # Determine alleles/ploidy
        if allele_idxs is not None and ploidy is not None:
            # Both specified: ploidy must match
            if len(allele_idxs) != ploidy:
                raise ValueError(
                    f"Ploidy of {ploidy} is not compatible with an "
                    f"allele_idxs list of {len(allele_idxs)} values"
                )
        elif allele_idxs is None:
            # allele_idxs not specified, use given ploidy or a default value of 2 to make a missing GT
            if ploidy is None:
                ploidy = 2
            allele_idxs = [
                MISSING_IDX,
            ] * ploidy

        # Ensure allele_idxs is a sorted tuple
        allele_idxs = tuple(sorted(allele_idxs))

        self.variant = variant
        self.allele_idxs = allele_idxs

        # Validate parameters
        for a in self.allele_idxs:
            if not self.variant.is_valid_allele_idx(a):
                raise ValueError(f"Invalid allele index for {self.variant}: {a}")

    @property
    def ploidy(self):
        """Number of sets of chromosomes"""
        return len(self.allele_idxs)

    def __str__(self):
        if all([a == MISSING_IDX for a in self.allele_idxs]):
            return "<Missing>"
        else:
            return "/".join([self.variant.alleles[a] for a in self.allele_idxs])

    def __repr__(self):
        return f"Genotype(variant={self.variant})[{str(self)}]"

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
