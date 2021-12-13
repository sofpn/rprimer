## Import data to test on
infile <- system.file("extdata", "example_alignment.txt", package = "rprimer")
testdata <- Biostrings::readDNAMultipleAlignment(infile)
testmat <- .consensusMatrix(testdata)

# consensusProfile =============================================================

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
    Biostrings::rowmask(testdata, invert = TRUE) <- 3
    prof <- consensusProfile(testdata)
    expect_s4_class(prof, "RprimerProfile")
    expect_true(all(prof$gaps == 0 | prof$gaps == 1))
    expect_true(all(prof$coverage == 1 | is.na(prof$coverage)))
})

test_that("consensusProfile works with a colmask", {
    Biostrings::colmask(testdata, invert = TRUE) <- 500:550
    prof <- consensusProfile(testdata)
    expect_s4_class(prof, "RprimerProfile")
})

# .consensusMatrix =============================================================

test_that(".consensusMatrix works", {
    expect_true(is.matrix(testmat))
    expect_equal(nrow(testmat), 6)
    expect_equal(rownames(testmat), c("A", "C", "G", "T", "-", "other"))
    expect_false(any(is.na(testmat)))
})

# .findMostCommonBase ==========================================================

test_that(".findMostCommonbase works", {
    testvec <- testmat[, 2]
    expect_equal(.findMostCommonBase(testvec), "-")
})

# .majorityConsensus ===========================================================

test_that(".majorityConsensus works", {
    testmat <- testmat[, 1:5]
    expect_equal(.majorityConsensus(testmat), rep("-", 5))
})

# .asIUPAC =====================================================================

test_that(".asIUPAC works", {
    expect_true(is.na(.asIUPAC("GC")))
    expect_equal(.asIUPAC("C,G"), "S")
    expect_equal(.asIUPAC("C,G,G"), "S")
    expect_equal(.asIUPAC("-"), "-")
    expect_true(is.na(.asIUPAC("X")))
})

# .dnaBasesOnly ================================================================

test_that(".dnaBasesOnly works", {
    dnaBases <- .dnaBasesOnly(testmat)[, 1:10]
    expect_equal(rownames(dnaBases), c("A", "C", "G", "T"))
    expect_true(
        all(colSums(dnaBases) > 0.9999999 & colSums(dnaBases < 1.00000001))
    )
})

# .iupacConsensus ==============================================================

test_that(".iupacConsensus works", {
    selection <- testmat[, 200:220]
    x <- selection[, 1:4]
    bases <- c("A", "C", "G", "T", "-")
    s <- x[rownames(x) %in% bases, , drop = FALSE]
    s <- apply(s, 2, function(x) x / sum(x))
    expect_equal(colSums(s), rep(1, 4), ignore_attr = TRUE)
})

# .nucleotideIdentity ==========================================================

test_that(".nucleotideIdentity works", {
    selection <- testmat[, 200:220]
    x <- selection[, 1:4]
    bases <- c("A", "C", "G", "T", "-")
    s <- x[rownames(x) %in% bases, , drop = FALSE]
    s <- apply(s, 2, function(x) x / sum(x))
    expect_equal(
        .nucleotideIdentity(x),
        c(max(s[, 1]), max(s[, 2]), max(s[, 3]), max(s[, 4]))
    )
})

# .coverage ====================================================================
