from typing import Optional, List


class Variant:
    """
    Information about a variant.

    Parameters
    ----------
    chromosome: str, optional
        None by default, through this usually is not desired.
    coordinate: int, optional
        (1-based, the default is 0 for none/unknown)
    variant_id: str, optional
        None by default
    alleles: List[str], optional
        List of possible alleles, empty by default

    Examples
    --------
    >>> variant = Variant('12', 112161652, 'rs12462', alleles=['C', 'T'])
    >>> print(variant)
    rs12462[chr=12;pos=112161652;2 alleles]
    """
    def __init__(self,
                 chromosome: Optional[str] = None,
                 coordinate: int = 0,
                 variant_id: Optional[str] = None,
                 alleles: Optional[List[str]] = None):
        self.chromosome = chromosome
        self.coordinate = coordinate
        self.variant_id = variant_id
        # Start with standard alleles - any additional (ins/del) are added to the end for sorting purposes
        self.alleles = ['A', 'C', 'G', 'T']
        if alleles is not None:
            self.alleles += [a for a in alleles if a not in self.alleles]

        # Validate the passed parameters
        if self.chromosome is not None and (';' in self.chromosome or ',' in self.chromosome):
            raise ValueError(f"The chromosome cannot contain ';' or ',': '{self.chromosome}'")
        if self.coordinate > ((2 ** 31) - 2):
            raise ValueError(f"The coordinate value may not exceed 2^31-2, {self.coordinate:,} was specified")
        if self.variant_id is not None and (';' in self.variant_id or ',' in self.variant_id):
            raise ValueError(f"The variant_id cannot contain ';' or ',': '{self.variant_id}'")
        if len(self.alleles) > 255:
            raise ValueError(f"{len(self.alleles):,} alleles were provided, the maximum supported number is 255.")

    def __str__(self):
        return f"{self.variant_id}[chr={self.chromosome};pos={self.coordinate};{len(self.alleles)} alleles]"

    def __repr__(self):
        return f"Variant(chromosome={self.chromosome}, coordinate={self.coordinate}," \
               f"variant_id={self.variant_id}, alleles={self.alleles})"

    def __eq__(self, other):
        if other.__class__ is not self.__class__:
            return NotImplemented
        return (self.chromosome == other.chromosome) &\
               (self.coordinate == other.coordinate) &\
               (self.variant_id == other.variant_id) & \
               (self.alleles == other.alleles)

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
        if len(self.alleles) < 255:
            self.alleles.append(allele)
        else:
            raise ValueError(f"Couldn't add new allele to {self}, 255 alleles max.")
        print(self.alleles)

    def get_allele_idx(self, allele: Optional[str], add: bool = False) -> int:
        """
        Get the integer value for an allele based on this variant's list of alleles,
        optionally adding it as a new allele

        Parameters
        ----------
        allele: str
        add: bool
            If `add` is True and the allele doesn't exist in the variant, add it
            If 'add' is False and the allele doesn't exist in the variant, return 255

        Returns
        -------
        int
            The integer value for the allele in the variant

        """
        if allele is None:
            return 255
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
                    raise ValueError(f"'{allele}' is not an allele in {self} a")
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
            True if the index value is valid (including 255 for None/Missing)
            False if the index value is not valid (negative or a nonexistent allele index)

        """
        if idx == 255:
            # None/Missing
            return True
        elif idx < 0:
            return False
        elif idx > (len(self.alleles)-1):
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
            True if the variant_id, chromosome, and coordinate all match.
            False otherwise.

        """
        if isinstance(other, Variant):
            return (self.variant_id == other.variant_id and
                    self.chromosome == other.chromosome and
                    self.coordinate == other.coordinate)
        else:
            return False

    def make_genotype(self, allele1: Optional[str] = None, allele2: Optional[str] = None) -> 'Genotype':
        """
        Create a Genotype object associated with this variant using the specified allele(s)

        Parameters
        ----------
        allele1: Optional[str]
        allele2: Optional[str]

        Returns
        -------
        Genotype
            A Genotype based on this variant with the specified alleles
        """
        a1 = self.get_allele_idx(allele1, add=True)
        a2 = self.get_allele_idx(allele2, add=True)
        return Genotype(self, a1, a2)

    def make_genotype_from_str(self, gt_str: str, sep: str = "/") -> 'Genotype':
        """
        Create a Genotype object associated with this variant using the specified allele string

        Parameters
        ----------
        gt_str: str
            A string containing two alleles separated by a delimiter character
        sep: str
            The delimiter used to separate alleles in the gt_str ('/' by default)

        Returns
        -------
        Genotype
            A Genotype based on this variant with the specified alleles
        """
        # Process Allele String
        alleles = gt_str.split(sep)
        if len(alleles) == 0:
            allele1 = None
            allele2 = None
        elif len(alleles) == 1:
            allele1 = alleles[0]
            allele2 = None
        elif len(alleles) == 2:
            allele1 = alleles[0]
            allele2 = alleles[1]
        else:
            raise ValueError("Can't process more that two alleles for a genotype")

        a1 = self.get_allele_idx(allele1, add=True)
        a2 = self.get_allele_idx(allele2, add=True)
        return Genotype(self, a1, a2)

    def make_genotype_from_plink_bits(self, plink_bits: str) -> 'Genotype':
        """
        Create a genotype from PLINK Bed file bits

        Parameters
        ----------
        plink_bits: str
            A string with allele indices as encoded in plink format, one of {'00', '01', '10', '11'}

        Returns
        -------
        Genotype
            A Genotype based on this variant with the specified alleles
        """
        # Process Allele String
        if plink_bits == '00':
            a1 = 0
            a2 = 0
        elif plink_bits == '01':
            a1 = 255
            a2 = 255
        elif plink_bits == '10':
            a1 = 0
            a2 = 1
        elif plink_bits == '11':
            a1 = 1
            a2 = 1
        else:
            raise ValueError(f"Invalid plink_bits: '{plink_bits}'")

        return Genotype(self, a1, a2)


class Genotype:
    """
    Genotype information associated with a specific variant.
    Defaults to using an anonymous variant with unknown alleles.
    Usually created with methods on ~Variant

    Parameters
    ----------
    variant: pandas_genomics.scalars.variant.Variant
    allele1: int
        The first allele encoded as an index into the variant allele list
    allele2: int
        The second allele encoded as an index into the variant allele list

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
    def __init__(self,
                 variant: Variant,
                 allele1: int = 255,
                 allele2: int = 255):
        self.variant = variant
        self.allele1 = allele1
        self.allele2 = allele2

        # Sort allele1 and allele2
        if self.allele1 > self.allele2:
            a1, a2 = self.allele2, self.allele1
            self.__setattr__('allele1', a1)
            self.__setattr__('allele2', a2)
        # Validate parameters
        if not self.variant.is_valid_allele_idx(self.allele1):
            raise ValueError(f"Invalid allele1 for {self.variant}: {self.allele1}")
        if not self.variant.is_valid_allele_idx(self.allele2):
            raise ValueError(f"Invalid allele2 for {self.variant}: {self.allele2}")

    def __str__(self):
        if self.allele1 == 255 and self.allele2 == 255:
            return "<Missing>"
        elif self.allele1 != 255 and self.allele2 == 255:
            return self.variant.alleles[self.allele1]
        elif self.allele1 != 255 and self.allele2 != 255:
            return f"{self.variant.alleles[self.allele1]}/{self.variant.alleles[self.allele2]}"

    def __repr__(self):
        return f"Genotype(variant={self.variant}, allele1={self.allele1}, allele2={self.allele2})"

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, other):
        if other.__class__ is not self.__class__:
            return NotImplemented
        return (self.variant == other.variant) & (self.allele1 == other.allele1) & (self.allele2 == other.allele2)

    def __lt__(self, other):
        if self.variant != other.variant:
            raise NotImplementedError("Can't compare different variants")
        else:
            # Compare allele index values for sorting
            # allele_1 is always <= allele_2 within a genotype
            a1_lt = self.allele1 < other.allele1
            a1_eq = self.allele1 == other.allele1
            a2_lt = self.allele2 < other.allele2
            return a1_lt | (a1_eq & a2_lt)

    def __gt__(self, other):
        if self.variant != other.variant:
            raise NotImplementedError("Can't compare different variants")
        else:
            # Compare allele index values for sorting
            # allele_1 is always <= allele_2 within a genotype
            a1_gt = self.allele1 > other.allele1
            a1_eq = self.allele1 == other.allele1
            a2_gt = self.allele2 > other.allele2
            return a1_gt | (a1_eq & a2_gt)

    def __le__(self, other):
        return (self < other) | (self == other)

    def __ge__(self, other):
        return (self > other) | (self == other)

    def is_missing(self) -> bool:
        """
        Returns
        -------
        bool
            True if the variant is missing (both alleles are None), otherwise False
        """
        return (self.allele1 == 255) and (self.allele2 == 255)
