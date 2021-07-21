from pathlib import Path

from pandas_genomics import io
from pandas_genomics.scalars import Region

DATA_DIR = Path(__file__).parent.parent / "data" / "bed"


def test_bed():
    input = DATA_DIR / "bed_test.bed"
    result = io.from_bed(input)
    assert list(result) == [
        Region("chr7", 127471197, 127472364, "Pos1"),
        Region("chr7", 127472364, 127473531, "Pos2"),
        Region("chr7", 127473531, 127474698, "Pos3"),
    ]
