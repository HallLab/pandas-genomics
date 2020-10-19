from pandas_genomics.scalars import Variant


def test_create_variant():
    variant = Variant("12", 112161652, "rs12462", ref="C", alt=["T"])
    print(variant)
