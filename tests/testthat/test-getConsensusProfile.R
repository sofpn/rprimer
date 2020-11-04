infile <- system.file("extdata", "example_alignment.txt", package = "rprimer")
testdata <- Biostrings::readDNAMultipleAlignment(infile)

test_that("getConsensusProfile returns an error when it should", {
    expect_error(getConsensusProfile(unclass(testdata)))
    expect_error(getConsensusProfile(testdata, iupacThreshold = FALSE))
    expect_error(getConsensusProfile(testdata, iupacThreshold = -0.1))
    expect_error(getConsensusProfile(testdata, iuacThreshold = 0.21))
})

test_that("getConsensusProfile works", {
    expect_s4_class(getConsensusProfile(testdata), "RprimerProfile")
})

testmatr <- .getConsensusMatrix(testdata)

test_that(".getConsensusMatrix works", {
    expect_equal(rownames(testmatr), c("A", "C", "G", "T", "-", "other"))
    expect_equal(min(colSums(testmatr)), 1)
    expect_equal(max(colSums(testmatr)), 1)
    expect_false(any(is.na(testmatr)))
})

test_that(".majorityConsensus works", {
    expect_equal(.majorityConsensus(testmatr[, 200:203]), rep("T", 4))
})

test_that(".asIUPAC works", {
    expect_true(is.na(.asIUPAC("GC")))
    expect_equal(.asIUPAC("C,G"), "S")
    expect_equal(.asIUPAC("C,G,G"), "S")
    expect_equal(.asIUPAC("C , G"), "S")
    expect_equal(.asIUPAC("-"), "-")
    expect_true(is.na(.asIUPAC("X")))
})

test_that(".iupacConsensus returns a warning when it should", {
    expect_warning(.iupacConsensus(testmatr[1:4, 1:10]))
})

test_that(".iupacConsensus works", {
    x <- testmatr[, 104:105]
    expect_equal(.iupacConsensus(x), c("H", "C"))
    expect_equal(.iupacConsensus(x, iupacThreshold = 0.03), c("Y", "C"))
})

test_that(".nucleotideIdentity works", {
    x <- testmatr[, 1:5]
    expect_equal(.nucleotideIdentity(x), rep(1, 5))
})

test_that(".shannonEntropy works", {
    x <- testmatr[, 1:5]
    expect_equal(.shannonEntropy(x), rep(0, 5))
})
