def test_loaded(vcf_test):
    """Validate the dataset"""
    # TODO: Add more assertions
    assert vcf_test.shape == (189, 115)
    var0 = vcf_test["Variant_0"].array.variant
    assert var0.ref == "CCCTAA"
    assert var0.score == 92
