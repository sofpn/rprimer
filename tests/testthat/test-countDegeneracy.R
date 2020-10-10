test_that(".countDegeneracy works", {
    expect_equal(.countDegeneracy("N"), 4)
    expect_equal(.countDegeneracy("NRC"), 4 * 2 * 1)
    expect_equal(.countDegeneracy("nrc-"), 4 * 2 * 1)
})
