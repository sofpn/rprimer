## Import data to test on
data("exampleRprimerOligo")
x <- examplleRprimerOligo

# assays =======================================================================

test_that("assays returns an error when it sould", {
    expect_error(assays(unclass(x)))
    expect_error(assays(x, lengthRange = c(39, 120)))
    expect_error(assays(x, lengthRange = c(40, 50001)))
    expect_error(assays(x, tmDiffPrimers = -0.1))
    expect_error(assays(x, tmDiffPrimers = 21))
    expect_error(assays(x, tmDiffPrimersProbe = -21))
    expect_error(assays(x, tmDiffPrimersProbe = 21))
})

test_that("assays works", {
    test <- assays(x)
    expect_s4_class(test, "RprimerAssay")
})

# .pairPrimers =================================================================

#test_that(".pairPrimers works", {
    #' .pairPrimers(x)
#})

# .combinePrimers ==============================================================

#test_that(".combinePrimers works", {

#})
