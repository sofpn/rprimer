test_that("getAlignmentProperties works", {
    data("exampleRprimerProfile")
    toTest <- exampleRprimerProfile
    expect_equal(length(majorityConsensus(toTest)), ncol(toTest))
    expect_equal(length(majorityConsensus(toTest[, 1:4])), 4)
    expect_true(is.character(iupacConsensus(toTest)))
    expect_false(any(grepl("^ACGT-", majorityConsensus(toTest))))
    expect_warning(iupacConsensus(toTest[1:2, ], threshold = 0.2))
    gaps <- gapFrequency(toTest)
    expect_true(is.double(gaps))
    expect_false(any(gaps < 0))
    expect_false(any(gaps > 1))
    identity <- nucleotideIdentity(toTest)
    expect_true(is.double(identity))
    expect_false(any(identity < 0))
    expect_false(any(identity > 1))
    expect_false(any(is.na(identity)))
    bases <- c("A", "C", "G", "T")
    s <- toTest[which(rownames(toTest) %in% bases), ]
    expect_equal(nucleotideIdentity(s), identity)
    entropy <- shannonEntropy(toTest)
    expect_true(is.double(entropy))
    expect_false(any(entropy < 0))
    expect_false(any(is.na(entropy)))
    expect_equal(shannonEntropy(s), entropy)
    properties <- getAlignmentProperties(toTest)
})

test_that("getAlignmentProperties returns an error when it should", {
    data("exampleRprimerProfile")
    toTest <- exampleRprimerProfile
    expect_error(iupacConsensus(toTest, threshold = 3))
    expect_error(iupacConsensus(toTest, threshold = NA))
    expect_error(iupacConsensus(unclass(toTest)))
    expect_error(gapFrequency(unclass(toTest)))
    expect_error(nucleotideIdentity(unclass(toTest)))
    expect_error(shannonEntropy(unclass(toTest)))
    expect_error(getAlignmentProperties(unclass(toTest)))
    expect_error(getAlignmentProperties(toTest, iupacThreshold = 3))
    expect_error(getAlignmentProperties(toTest, iupacThreshold = NA))
})

test_that("asIUPAC returns an error when it should", {
    expect_error(asIUPAC(c("C,G", "G,T")))
})

test_that("asIUPAC works", {
    expect_true(is.na(asIUPAC("GC")))
    expect_equal(asIUPAC("C,G"), "S")
    expect_equal(asIUPAC("C,G,G"), "S")
    expect_equal(asIUPAC("C , G"), "S")
    expect_equal(asIUPAC("-"), "-")
})
