from pathlib import Path
from typing import Union

import pandas as pd

from pandas_genomics.scalars import Region


def from_bed(filename: Union[str, Path]):
    """
    Yields genomic regions from a bed file as Region scalars

    Parameters
    ----------
    filename: str or Path
        bed file

    Returns
    -------
    List[Region]
    """
    with open(filename, "r") as input:
        for idx, line in enumerate(input):
            if (
                line.startswith("browser")
                or line.startswith("track")
                or line.startswith("#")
            ):
                continue
            line = line.split("\t")
            if len(line) < 3:
                raise ValueError(
                    f"Expected at least 3 tab-delimited fields, found {len(line)} in '{line}'."
                )
            # Note that positions need to be incremented by 1 to go from 0-based to 1-based.
            yield Region(line[0], int(line[1]) + 1, int(line[2]) + 1, line[3])
