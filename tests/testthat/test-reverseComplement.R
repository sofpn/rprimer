test_that(".reverseComplement works", {
    expect_equal(.reverseComplement("acg"), "CGT")
    expect_equal(.reverseComplement("ACG"), "CGT")
    expect_equal(.reverseComplement("yrn"), "NYR")
    expect_equal(.reverseComplement("-a"), "T-")
})

test_that(".reverseComplement returns an error when it should", {
    expect_error(.reverseComplement(2))
    expect_error(.reverseComplement(c("g", "t")))
    expect_error(.reverseComplement("gxttn"))
})
