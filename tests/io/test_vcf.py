def test_loaded(vcf_test):
    """Validate the dataset"""
    # TODO: Add more assertions
    assert vcf_test.shape == (189, 115)
