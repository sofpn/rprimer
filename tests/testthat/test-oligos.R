## Import data to test on
data("exampleRprimerProfile")
x <- exampleRprimerProfile

# oligos =======================================================================

test_that("oligos returns an error when it should", {
    expect_error(oligos(unclass(x)))
    expect_error(oligos(x, maxGapFrequency = -1))
    expect_error(oligos(x, maxGapFrequency = 1.1))
    expect_error(oligos(x, maxGapFrequency = -1))
    expect_error(oligos(x, lengthPrimer = 14))
    expect_error(oligos(x, lengthPrimer = 41))
    expect_error(oligos(x, maxDegeneracyPrimer = 65))
    expect_error(oligos(x, maxDegeneracyPrimer = 0))
    expect_error(oligos(x, gcClampPrimer = 1))
    expect_error(oligos(x, avoidThreeEndRunsPrimer = "true"))
    expect_error(oligos(x, gcPrimer = c(-0.1, 1)))
    expect_error(oligos(x, gcPrimer = c(0, 1.1)))
    expect_error(oligos(x, tmPrimer = c(19, 90)))
    expect_error(oligos(x, tmPrimer = c(20, 91)))
    expect_error(oligos(x, concPrimer = 19))
    expect_error(oligos(x, concPrimer = 2001))
    expect_error(oligos(x, designStrategyPrimer = ""))
    expect_error(oligos(x, probe = 1))
    expect_error(oligos(x, lengthProbe = 14))
    expect_error(oligos(x, lengthProbe = 41))
    expect_error(oligos(x, maxDegeneracyProbe = 0))
    expect_error(oligos(x, maxDegeneracyProbe = 65))
    expect_error(oligos(x, avoidFiveEndGProbe = 1))
    expect_error(oligos(x, gcProbe = c(-0.1, 1)))
    expect_error(oligos(x, gcProbe = c(0, 1.1)))
    expect_error(oligos(x, tmProbe = c(19, 90)))
    expect_error(oligos(x, tmProbe = c(20, 91)))
    expect_error(oligos(x, concProbe = 19))
    expect_error(oligos(x, concProbe = 2001))
    expect_error(oligos(x, concNa = 0.009))
    expect_error(oligos(x, concNa = 1.1))
})

test_that("oligos works", {
    test <- oligos(x)
    expect_s4_class(test, "RprimerOligo")
    test <- oligos(x, probe = FALSE)
    expect_true(all(test$type == "primer"))
    test <- oligos(
        x,
        designStrategyPrimer = "mixed",
        maxDegeneracyPrimer = 1,
        maxDegeneracyProbe = 2,
        lengthPrimer = 25:30,
        lengthProbe = 18
    )
    primers <- test[test$type == "primer", ]
    expect_true(all(primers$method[primers$fwd] == "mixedFwd"))
    expect_true(all(primers$method[primers$rev] == "mixedRev"))
    expect_true(max(primers$degeneracy) == 1)
    expect_true(all(primers$length >= 25 & primers$length <= 30))
    probes <- test[test$type == "probe", ]
    expect_true(all(probes$method == "ambiguous"))
    expect_true(all(probes$length == 18))
    expect_true(max(probes$degeneracy) == 2)
    expect_error(oligos(x[1:100, ]))
    expect_error(oligos(x[1:1000, ], maxDegeneracyProbe = 1))

    ## test with only one sequence
    infile <- system.file("extdata", "example_alignment.txt", package = "rprimer")
    testdata <- Biostrings::readDNAMultipleAlignment(infile)
    Biostrings::rowmask(testdata, invert = TRUE) <- 3
    prof <- consensusProfile(testdata)
    test <- oligos(prof[1:100, ])
    expect_s4_class(test, "RprimerOligo")
})

# .nmers =======================================================================

test_that(".nmers works", {
    seq <- sample(c("A", "C", "G", "T"), 100, replace = TRUE)
    nmer <- .nmers(seq, n = 4)
    expect_equal(ncol(nmer), 4)
    expect_equal(seq[1:4], nmer[1, ])
    expect_equal(seq[2:5], nmer[2, ])
    expect_equal(seq[(length(seq) - 3):length(seq)], nmer[nrow(nmer), ])

    ## Confirm that it returns a matrix even if it is only one row
    nmer <- .nmers(c("A", "C", "G"), 3)
    expect_true(is.matrix(nmer))
    expect_equal(dim(nmer), c(1, 3))
})

# .countDegeneracy =============================================================

test_that(".countDegeneracy works", {
    seq <- c("A", "C", "R", "N", "Y")
    degen <- .countDegeneracy(seq)
    table <- rprimer:::lookup$degeneracy
    expect_equal(degen, unname(table["R"] * table["N"] * table["Y"]))
    expect_equal(.countDegeneracy(c("A", "C", "G", "T")), 1)
})

# .splitAndPaste ===============================================================

test_that(".splitAndPaste works", {
    first <- t(matrix(rep(0, 12)))
    second <- t(matrix(rep(1, 12)))
    test <- .splitAndPaste(first, second)
    testRev <- .splitAndPaste(first, second, rev = TRUE)

    expect_true(is.matrix(test))
    expect_equal(dim(test), c(1, 12))
    expect_equal(dim(testRev), c(1, 12))

    ## If fwd, first (0) should be two thirds (= 8)
    ## and second (1) one third (= 4)
    expect_equal(sum(test[, 1:8]), 0)
    expect_equal(sum(test[9:12]), 4)

    ## If rev, first (0) should be one third (= 4)
    ## and second (1) two thirds (= 8)
    expect_equal(sum(testRev[, 1:4]), 0)
    expect_equal(sum(testRev[5:12]), 8)

    # Confirm that it works for numbers not dividable with three
    first <- matrix(rep(0, 17 * 4), ncol = 17, nrow = 4)
    second <- matrix(rep(1, 17 * 4), ncol = 17, nrow = 4)
    test <- .splitAndPaste(first, second)
    testRev <- .splitAndPaste(first, second, rev = TRUE)
    expect_equal(
        ncol(test[, colSums(test) == 0]),
        ncol(testRev[, colSums(testRev) == nrow(testRev)])
    )

    testSeparate <- .splitAndPaste(first, second, combine = FALSE)
    expect_true(is.list(testSeparate))
    expect_true(is.matrix(testSeparate[[1]]))
    expect_true(is.matrix(testSeparate[[2]]))
})

# .generateAmbiguousOligos =====================================================

test_that(".generateAmbiguousOligos works", {
    test <- .generateAmbiguousOligos(x, lengthOligo = 18)
    expect_equal(ncol(test$iupacSequence), 18)
    expect_equal(unique(test$length), 18)
    expect_equal(unique(test$roiStart), min(x$position))
    expect_equal(unique(test$roiEnd), max(x$position))
    expect_equal(
        x$iupac[test$start[4]:test$end[4]], test$iupacSequence[4, ]
    )
    expect_equal(
        test$degeneracy[4],
        .countDegeneracy(x$iupac[test$start[4]:test$end[4]])
    )
    expect_equal(
        test$gapFrequency[4], max(x$gaps[test$start[4]:test$end[4]])
    )
    expect_equal(
        test$coverage[4], mean(x$coverage[test$start[4]:test$end[4]])
    )
})

# .mergeLists ==================================================================

test_that(".mergeLists works", {
    l1 <- list("M" = matrix(rep(1, 10)), "V" = rep("A", 10))
    l2 <- list("M" = matrix(rep(2, 10)), "V" = rep("B", 10))
    l <- .mergeLists(l1, l2)
    expect_equal(names(l), names(l1))
    expect_true(is.matrix(l[[1]]))
    expect_true(is.vector(l[[2]]))
})

# .filterOligos ================================================================

test_that(".filterOligos works", {
    x <- .generateAmbiguousOligos(x)
    test <- .filterOligos(x, maxDegeneracy = 32, maxGapFrequency = 0)
    expect_true(all(test$degeneracy <= 32))
    expect_true(all(test$gapFrequency == 0))
    expect_true(all(!is.na(test$majoritySequence)))
    expect_true(all(test$majoritySequence != "-"))
})

# .expandDegenerates ===========================================================

test_that(".expandDegenerates works", {
    seq <- c("A", "R", "T", "T", "N", "G")
    degen <- .expandDegenerates(seq)
    nDegen <- .countDegeneracy(seq)
    expect_equal(nrow(degen), nDegen)
    seq2 <- c("A", "C", "G", "T")
    degen2 <- .expandDegenerates(seq2)
    expect_true(is.matrix(degen2))
    expect_equal(nrow(degen2), 1)
})

# .oligoMatrix =================================================================

test_that(".oligoMatrix works", {
    y <- .generateAmbiguousOligos(x[5000:6000, ])
    y <- .filterOligos(y, maxDegeneracy = 8)
    oligoList <- apply(y$iupacSequence, 1, .expandDegenerates)
    expect_true(all(vapply(oligoList, is.matrix, logical(1))))
    expect_equal(
        ncol(y$iupacSequence), unique(vapply(oligoList, ncol, integer(1)))
    )
    oligoMatrix <- .oligoMatrix(oligoList)
    expect_equal(length(unique(rownames(oligoMatrix))), length(oligoList))

    ## Test how it works when the maximum degeneracy is one
    y <- .generateAmbiguousOligos(x[5000:6000, ])
    y <- .filterOligos(y, maxDegeneracy = 1)
    oligoList <- apply(y$iupacSequence, 1, .expandDegenerates)
    if (!is.list(oligoList)) {
        oligoList <- t(oligoList)
        oligoList <- lapply(seq_len(nrow(oligoList)), function(i) {
            oligoList[i, , drop = FALSE]
        })
    }
    expect_equal(
        ncol(y$iupacSequence), unique(vapply(oligoList, ncol, integer(1)))
    )
    oligoMatrix <- .oligoMatrix(oligoList)
    expect_equal(length(unique(rownames(oligoMatrix))), length(oligoList))

    ## Test if there is only one oligo as input
    y <- lapply(y, function(x) {
        if (is.matrix(x)) x[1, , drop = FALSE] else x[1]
    })
    oligoList <- apply(y$iupacSequence, 1, .expandDegenerates)
    if (!is.list(oligoList)) oligoList <- list(t(oligoList))
    oligoMatrix <- .oligoMatrix(oligoList)
    expect_equal(length(unique(rownames(oligoMatrix))), length(oligoList))

    ## Test if degeneracy is 64 (max)
    y <- .generateAmbiguousOligos(x[5000:6000, ])
    y <- .filterOligos(y, maxDegeneracy = 64)
    oligoList <- apply(y$iupacSequence, 1, .expandDegenerates)
    expect_true(all(vapply(oligoList, is.matrix, logical(1))))
    expect_equal(
        ncol(y$iupacSequence), unique(vapply(oligoList, ncol, integer(1)))
    )
    oligoMatrix <- .oligoMatrix(oligoList)
    expect_equal(length(unique(rownames(oligoMatrix))), length(oligoList))
})

# .reverseComplement ===========================================================

test_that(".reverseComplement works", {
    seq <- t(matrix(c("A", "C", "G")))
    rc <- .reverseComplement(seq)
    expect_true(is.matrix(rc))
    expect_equal(dim(rc), c(1, 3))
    expect_equal(rc, c("C", "G", "T"), ignore_attr = TRUE)
})

# .gcClamp =====================================================================

test_that(".gcClamp works", {
    seq <- c("A", "C", "G", "T", "G", "C", "T", "A")
    gc <- as.integer(seq == "C" | seq == "G")
    gc <- t(matrix(gc))
    expect_true(.gcClamp(gc))
    seq <- c("A", "C", "G", "T", "T", "C", "T", "A")
    gc <- ifelse(seq == "C" | seq == "G", 1, 0)
    gc <- t(matrix(gc))
    expect_false(.gcClamp(gc))
    seq <- c("A", "C", "G", "T", "C", "C", "G", "C")
    gc <- ifelse(seq == "C" | seq == "G", 1, 0)
    gc <- t(matrix(gc))
    expect_false(.gcClamp(gc))
    seq <- t(matrix(c("A", "C", "G", "T", "C", "C", "G", "C")))
    rc <- .reverseComplement(seq)
    gcRc <- ifelse(rc == "C" | rc == "G", 1, 0)
    gcRc <- t(matrix(gcRc))
    expect_true(.gcClamp(gc, rev = TRUE))
    expect_true(.gcClamp(gcRc))
    gc <- matrix(rep(c(1, 0, 1, 1, 0, 0, 0, 0, 0), 8), ncol = 9, byrow = TRUE)
    expect_equal(sum(.gcClamp(gc)), 0)
    expect_equal(sum(.gcClamp(gc, rev = TRUE)), nrow(gc))
})

# .endRuns =====================================================================

test_that(".endRuns works", {
    seq <- t(matrix(c("A", "T", "C", "C", "C")))
    expect_true(.endRuns(seq))
    seq <- t(matrix(c("A", "T", "T", "C", "C")))
    expect_false(.endRuns(seq))
    expect_false(.endRuns(seq, rev = TRUE))
    seq <- t(matrix(c("G", "G", "G", "T", "A")))
    expect_true(.endRuns(seq, rev = TRUE))
    seq <- matrix(rep(c("G", "G", "G", "T", "A"), 10), ncol = 5, byrow = TRUE)
    expect_equal(sum(.endRuns(seq)), 0)
    expect_equal(sum(.endRuns(seq, rev = TRUE)), 10)
})

# .repeats =====================================================================

test_that(".repeats works", {
    mono <- "CTTTTTA"
    di <- "CTCTCTCTCTA"
    monodi <- "CTCTCTCTCTACTTTTTA"
    expect_true(.repeats(mono))
    expect_true(.repeats(di))
    expect_true(.repeats(monodi))
    expect_equal(sum(.repeats(c(di, mono))), 2)
})

# .allVariants ==============================================================

test_that(".allVariants works", {
    y <- .generateAmbiguousOligos(x[5000:6000, ])

    ## Make sure it works if max degeneracy is one
    y <- .filterOligos(y, maxDegeneracy = 1)
    all <- .allVariants(y)
    expect_true(is.list(all))
    expect_equal(
        all$sequence[[1]], paste(y$iupacSequence[1, ], collapse = "")
    )
    expect_equal(
        all$sequence[[length(all$sequence)]],
        paste(y$iupacSequence[nrow(y$iupacSequence), ], collapse = "")
    )
    expect_equal(unique(vapply(all$sequence, length, integer(1))), 1)

    ## If degeneracy is high
    y <- .generateAmbiguousOligos(x[100:150, ])
    y <- .filterOligos(y, maxDegeneracy = 16)
    all <- .allVariants(y)
    nVariants <- lapply(all, function(x) vapply(x, length, integer(1)))
    nVariants <- do.call("rbind", nVariants)
    expect_true(all(apply(nVariants, 2, function(x) length(unique(x))) == 1))
    expect_true(all(nVariants <= 16))
})

# .meanRange ===================================================================

test_that(".meanRange works", {
    x <- .allVariants(.filterOligos(.generateAmbiguousOligos(x)))
    test <- .meanRange(x)
    range <- function(x) max(x) - min(x)
    expect_equal(test$gcContentMean, vapply(x$gcContent, mean, double(1L)))
    expect_equal(test$tmPrimerMean, vapply(x$tmPrimer, mean, double(1L)))
    expect_equal(test$tmProbeMean, vapply(x$tmProbe, mean, double(1L)))
    expect_equal(test$gcContentRange, vapply(x$gcContent, range, double(1L)))
    expect_equal(test$tmPrimerRange, vapply(x$tmPrimer, range, double(1L)))
    expect_equal(test$tmProbeRange, vapply(x$tmProbe, range, double(1L)))
})

# .makeOligoDf =================================================================

test_that(".makeOligoDf works", {
    x <- .filterOligos(.generateMixedOligos(x))
    test <- .makeOligoDf(x)
    expect_true(is.data.frame(test))
})

# .designAmbiguousOligos =======================================================

test_that(".designAmbiguousOligos works", {
    test <- .designAmbiguousOligos(
        x,
        maxDegeneracyPrimer = 4,
        maxGapFrequency = 0,
        lengthPrimer = 15:18,
        probe = FALSE
    )
    expect_true(is.data.frame(test))
    expect_true(all(test$degeneracy <= 4))
    expect_true(all(test$method == "ambiguous"))
    expect_true(
        all(test$length >= 15 & test$length <= 18)
    )
    expect_error(.designAmbiguousOligos(x[1:50, ], maxDegeneracy = 1))
})

# .isWithinRange ===============================================================

test_that(".isWithinRange works", {
    x <- .allVariants(.filterOligos(.generateAmbiguousOligos(x)))
    x$gcContent <- x$gcContent[1]
    test <- unlist(.isWithinRange(x$gcContent, c(0.5, 0.6)))
    expect_equal(test, c(TRUE, TRUE, TRUE, FALSE))
})

# .isValid =====================================================================

test_that(".isValid works", {
    x <- .designAmbiguousOligos(x)
    gcInRange <- .isWithinRange(x$gcContent, c(0.4, 0.6))
    x <- cbind(x, data.frame(cbind(gcInRange)))
    check <- .convertToMatrices(x[c("gcInRange", "repeats", "threeEndRunsFwd")])
    test <- .isValid(check, rowThreshold = 1, colThreshold = 1)
    expect_equal(nrow(x), length(test))
})

# .filterPrimers ===============================================================

# .checkAllProbeVariants =======================================================

# .filterProbes ================================================================

# .beautifyOligos ==============================================================
