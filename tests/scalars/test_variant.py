from pandas_genomics import Variant


def test_create_variant():
    variant = Variant('12', 112161652, 'rs12462', alleles=['C', 'T'])
    print(variant)
    print(variant == variant)
