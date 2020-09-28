test_that("countDegeneracy works", {
    expect_equal(countDegeneracy("N"), 4)
    expect_equal(countDegeneracy("NRC"), 4 * 2 * 1)
    expect_equal(countDegeneracy("nrc-"), 4 * 2 * 1)
})

test_that("countDegeneracy returns an error when it should", {
    expect_error(countDegeneracy(c("A", "G")))
    expect_error(countDegeneracy("XGT"))
    expect_error(countDegeneracy(0))
})
