import pytest


def test_loaded_small(plink_small):
    """Validate the small dataset"""
    # TODO: Add more assertions
    assert plink_small.shape == (150, 3020)


@pytest.mark.slow
def test_loaded_medium(plink_medium):
    """Validate the medium dataset"""
    # TODO: Add more assertions
    assert plink_medium.shape == (600, 45100)
