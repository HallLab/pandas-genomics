import pytest

from pandas_genomics.scalars import Region, Variant


@pytest.mark.parametrize(
    "chromosome,start,end",
    [
        ("chr1", 1, 2),
        pytest.param(
            "chr1",
            0,
            1,
            marks=pytest.mark.xfail(raises=ValueError, strict=True, reason="Start < 1"),
        ),
        pytest.param(
            "chr1",
            5,
            5,
            marks=pytest.mark.xfail(
                raises=ValueError, strict=True, reason="Start == End"
            ),
        ),
        pytest.param(
            "chr1",
            6,
            5,
            marks=pytest.mark.xfail(
                raises=ValueError, strict=True, reason="Start > End"
            ),
        ),
        pytest.param(
            3,
            6,
            5,
            marks=pytest.mark.xfail(
                raises=TypeError, strict=True, reason="chromosome not string"
            ),
        ),
    ],
)
def test_create_region(chromosome, start, end):
    region = Region(chromosome, start, end)


@pytest.mark.parametrize(
    "variant,region,result",
    [
        (
            Variant(chromosome="chr1", position=1),
            Region(chromosome="chr1", start=1, end=2),
            True,
        ),
        (
            Variant(chromosome="chr1", position=1),
            Region(chromosome="chr2", start=1, end=2),
            False,
        ),
        (
            Variant(chromosome="chr1", position=1),
            Region(chromosome="1", start=1, end=2),
            False,
        ),
        (
            Variant(chromosome="chr1", position=99),
            Region(chromosome="chr1", start=98, end=99),
            False,
        ),
    ],
)
def test_contains_variant(variant, region, result):
    assert region.contains_variant(variant) == result
