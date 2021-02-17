## Import data to test on
data("exampleRprimerProfile")
data("exampleRprimerOligo")

# RprimerProfile ===============================================================

test_that("RprimerProfile validation works", {
    expect_error(RprimerProfile(10))

    ## Should be subsettable by row, but not by column
    expect_error(exampleRprimerProfile[, 1:2])
    expect_s4_class(exampleRprimerProfile[1:10, ], "RprimerProfile")
})

# RprimerOligo =================================================================

test_that("RprimerOligo validation works", {
    expect_error(RprimerOligo(10))

    ## Should be subsettable by row, but not by column
    expect_error(exampleRprimerOligo[, 1:2])
    expect_s4_class(exampleRprimerOligo[1:10, ], "RprimerOligo")
})

# RprimerAssay =================================================================

test_that("RprimerAssay validation works", {
    expect_error(RprimerAssay(10))

    ## Should be subsettable by row, but not by column
    expect_error(exampleRprimerAssay[, 1:2])
    expect_s4_class(exampleRprimerAssay[1:10, ], "RprimerAssay")
})
