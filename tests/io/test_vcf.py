def test_loaded(vcf_test):
    """Validate the dataset"""
    # TODO: Add more assertions
    assert vcf_test.shape == (189, 115)
    var0 = vcf_test.iloc[:, 0]
    assert var0.genomics.variant.ref == "CCCTAA"
    assert var0.genomics.variant.score == 92
    filt_v0 = var0.loc[var0.genomics.gt_scores > 10]
    assert len(filt_v0) == 23
