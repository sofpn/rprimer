test_that(".gcContent works", {
    expect_equal(.gcContent("ACTCTC"), 0.5)
    expect_equal(.gcContent("ACTCTC--"), 0.5)
    expect_equal(.gcContent("actctc--"), 0.5)
})

