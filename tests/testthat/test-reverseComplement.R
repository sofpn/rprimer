test_that(".reverseComplement works", {
    expect_equal(.reverseComplement("acg"), "CGT")
    expect_equal(.reverseComplement("ACG"), "CGT")
    expect_equal(.reverseComplement("yrn"), "NYR")
    expect_equal(.reverseComplement("-a"), "T-")
})


