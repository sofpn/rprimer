data("exampleRprimerOligo")
testdata <- exampleRprimerOligo

test_that("getAssays works", {
    assays <- getAssays(testdata)
    expect_s4_class(assays, "RprimerAssay")
    probeBegin <- assays$startPr - assays$endFwd
    probeEnd <- assays$startRev - assays$endPr
    expect_false(any(probeBegin < 2))
    expect_false(any(probeEnd < 2))
})

test_that("getAssays returns an error when it should", {
    expect_error(getAssays(unclass(exampleRprimerOligo)))
})

test_that(".gContent works", {
    expect_equal(.gContent("GGGGTTTT"), 0.5)
})

test_that(".findAllPrimerCombinations works", {
    out <- .findAllPrimerCombinations(as.data.frame(testdata))
    startValidFwd <- testdata$start[!is.na(testdata$majority)]
    startValidRev <- testdata$start[!is.na(testdata$majorityRc)]
    expect_equal(length(setdiff(out$startFwd, startValidFwd)), 0L)
    expect_equal(length(setdiff(out$startRev, startValidRev)), 0L)
    expect_false(any(is.na(out$majorityFwd)))
    expect_false(any(is.na(out$majorityRevRc)))
})

test_that(".combinePrimers works", {
    assays <- .combinePrimers(
        as.data.frame(testdata),
        length = 60:200, maxTmDifferencePrimers = 2
    )
    expect_true(all(assays$tmDifferencePrimer <= 2))
    expect_true(all(assays$ampliconLength >= 60))
    expect_true(all(assays$ampliconLength <= 200))
})

test_that(".combinePrimers returns an error when it should", {
    expect_error(
        .combinePrimers(as.data.frame(testdata), maxTmDifferencePrimers = -0.1)
    )
    expect_error(
        .combinePrimers(as.data.frame(testdata), maxTmDifferencePrimers = 31)
    )
    expect_error(.combinePrimers(as.data.frame(testdata), length = 39:100))
    expect_error(.combinePrimers(as.data.frame(testdata), length = 40:5001))
    expect_error(.combinePrimers(as.data.frame(testdata[1, ])))
})

# test_that(".getAllProbeCandidates work", {
#  assays <- .
# })


# test_that(".addProbes works", {
# })

test_that(".addProbes returns an error when it should", {
    assays <- .combinePrimers(as.data.frame(testdata))
    probes <- as.data.frame(testdata[testdata$type == "probe", ])
    expect_error(
        .addProbes(assays, probes, tmDifferencePrimersProbe = c(-21, 10))
    )
    expect_error(
        .addProbes(assays, probes, tmDifferencePrimersProbe = c(-20, 21))
    )
    expect_error(.addProbes(assays[1, ], probes[1, ]))
})
