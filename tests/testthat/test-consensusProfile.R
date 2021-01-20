infile <- system.file("extdata", "example_alignment.txt", package = "rprimer")
testdata <- Biostrings::readDNAMultipleAlignment(infile)

test_that("consensusProfile returns an error when it should", {
    expect_error(consensusProfile(unclass(testdata)))
    expect_error(consensusProfile(testdata, -0.1))
    expect_error(consensusProfile(testdata, 0.21))
})

test_that("consensusProfile works", {
    expect_s4_class(consensusProfile(testdata), "RprimerProfile")
})

test_that("consensusProfile works with a rowmask", {
    ## Select only one sequence
    Biostrings::rowmask(testdata, invert = TRUE) <- 1
    prof <- consensusProfile(testdata)
    expect_s4_class(prof, "RprimerProfile")
    expect_true(all(prof$identity == 1))
    expect_true(all(prof$entropy == 0))
    expect_true(all(prof$gaps == 0 | prof$gaps == 1))
})

test_that("consensusProfile works with a colmask", {
    Biostrings::colmask(testdata, invert = TRUE) <- 500:550
    prof <- consensusProfile(testdata)
    expect_s4_class(prof, "RprimerProfile")
    expect_equal(prof$position, 1:51)
})

testmat <- .consensusMatrix(testdata)

test_that(".consensusMatrix works", {
    expect_equal(rownames(testmat), c("A", "C", "G", "T", "-", "other"))
    expect_equal(nrow(testmat), 6)
    expect_true(all(colSums(testmat), 1))
    expect_false(any(is.na(testmat)))
})

selection <- testmat[, 200:220]

test_that(".majorityConsensus works", {
    expect_equal(.majorityConsensus(selection[, 1:4]), c("T", "A", "T", "T"))
})

test_that(".asIUPAC works", {
    expect_true(is.na(.asIUPAC("GC")))
    expect_equal(.asIUPAC("C,G"), "S")
    expect_equal(.asIUPAC("C,G,G"), "S")
    expect_equal(.asIUPAC("-"), "-")
    expect_true(is.na(.asIUPAC("X")))
})

test_that(".iupacConsensus returns a warning when it should", {
    expect_warning(.iupacConsensus(selection[2:3, ], ambiguityThreshold = 0.2))
})

x <- selection[, 1:4]
bases <- c("A", "C", "G", "T", "-")
s <- x[rownames(x) %in% bases, , drop = FALSE]
s <- apply(s, 2, function(x) x / sum(x))

test_that(".iupacConsensus works", {
    expect_equal(.iupacConsensus(x), c("N", "A", "T", "H"))
    expect_equal(.iupacConsensus(x, 0.05), c("Y", "A", "T", "H"))
    expect_equivalent(colSums(s), rep(1, 4))
})

test_that(".nucleotideIdentity works", {
    expect_equal(
        .nucleotideIdentity(x),
        c(max(s[, 1]), max(s[, 2]), max(s[, 3]), max(s[, 4]))
    )
})

test_that(".shannonEntropy works", {
    expect_equal(.shannonEntropy(x), .shannonEntropy(s))
    entropy <- 0.035*log2(0.035) + 0.240*log2(0.240) + 0.045*log2(0.045) + 0.680*log2(0.680)
    entropy <- abs(entropy)
    expect_equal(.shannonEntropy(s[, 1, drop = FALSE]), entropy)
})

test_that(".coverage works", {
    expect_equivalent(
        .coverage(x[, 1, drop = FALSE], 0.2), 1 - (0.035 + 0.045)
    )
})
