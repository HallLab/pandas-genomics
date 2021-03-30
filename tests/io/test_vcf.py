def test_loaded(vcf_test):
    """Validate the dataset"""
    # TODO: Add more assertions
    assert vcf_test.shape == (189, 115)
    var0 = vcf_test["Variant_0"]
    assert var0.genotype.variant.ref == "CCCTAA"
    assert var0.genotype.variant.score == 92
    filt_v0 = var0.loc[var0.genotype.gt_scores > 10]
    assert len(filt_v0) == 23
