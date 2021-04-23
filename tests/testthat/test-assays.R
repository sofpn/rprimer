## Import data to test on
data("exampleRprimerOligo")
x <- exampleRprimerOligo

# assays =======================================================================

test_that("assays returns an error when it sould", {
    expect_error(assays(unclass(x)))
    expect_error(assays(x, lengthRange = c(39, 120)))
    expect_error(assays(x, lengthRange = c(40, 50001)))
})

test_that("assays works", {
    test <- assays(x)
    expect_s4_class(test, "RprimerAssay")
})

# .pairPrimers =================================================================

# .combinePrimers ==============================================================

test_that(".combinePrimers works", {
    test <- .combinePrimers(x,
        lengthRange = c(100, 2000)
    )
    expect_true(all(test$ampliconLength >= 100))
    expect_true(all(test$ampliconLength <= 2000))
    expect_equal(test$end - test$start, test$ampliconLength - 1)
    expect_error(.combinePrimers(x[1, ]))
})

# .identifyProbes ==============================================================

test_that(".identifyProbes works", {
    assays <- .combinePrimers(x)
    test <- .identifyProbes(assays, x[x$type == "probe", ])
    expect_equal(nrow(assays), length(test))
})

# .extractProbes ===============================================================

test_that(".extractProbes works", {
    assays <- .combinePrimers(x)
    probes <- .identifyProbes(assays, x[x$type == "probe", ])
    test <- .extractProbes(assays, probes)
    expect_true(all(test$startPr - test$endFwd >= 1))
    expect_true(all(test$startRev - test$endPr >= 1))
})

# .addProbes ===================================================================

# .beautifyPrimers =============================================================

# .beautifyProbes ==============================================================
