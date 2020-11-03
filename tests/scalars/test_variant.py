from pandas_genomics.scalars import Variant


def test_create_variant():
    variant = Variant("12", 112161652, "rs12462", ref="C", alt=["T"])
    assert variant.alleles == ["C", "T"]


def test_methods():
    variant = Variant("12", 112161652, "rs12462", ref="C", alt=["T"])
    variant_also = Variant("12", 112161652, "rs12462", ref="C", alt=["T"])
    assert variant.is_same_variant(variant_also)
    # Get Allele Index
    assert variant.get_allele_idx("T") == 1
    assert variant.get_allele_idx("G", add=True) == 2
    assert len(variant.alleles) == 3
    # Add Allele
    variant.add_allele("GT")
    assert len(variant.alleles) == 4
    # Is Valid Allele Index
    assert variant.is_valid_allele_idx(1)
    assert not variant.is_valid_allele_idx(10)
    # Same variant despite adding additional alleles
    assert variant.is_same_variant(variant_also)
    # But variant not equal
    assert not variant == variant_also
