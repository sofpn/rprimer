# Skriv stort test för olika fall

data("exampleRprimerProfile")
x <- exampleRprimerProfile
y <- .generateOligos(x[5000:6000, ])


test_that("oligos works", {
    z <- oligos(x[5000:6000, ],
                probe = FALSE,
                lengthPrimer = 16:20,
                maxDegeneracyPrimer = 2,
                tmRangePrimer = c(50, 60),
                gcRangePrimer = c(0.45, 0.75))
    expect_true(all(z$type == "primer"))
    expect_true(all(unlist(z$gcContent) >= 0.45))
    expect_true(all(unlist(z$gcContent) <= 0.75))
    expect_true(all(z$degeneracy <= 2))
    di <- "(AT){4,}|(TA){4,}|(AC){4,}|(CA){4,}|(AG){4,}|(GA){4,}|(GT){4,}|(TG){4,}|(CG){4,}|(GC){4,}|(CT){4,}|(TC){4,}|)"
    mono <- "([A-Z])\\1\\1\\1\\1"
    expect_false(any(grepl(di, unlist(z$sequence))))
    expect_false(any(grepl(mono, unlist(z$sequence))))
    expect_true(all(unlist(z$tm) >= 50))
    expect_true(all(unlist(z$tm) <= 60))
    expect_equal(z$tmMean, vapply(z$tm, mean, double(1)))
    expect_equal(z$tmRange, vapply(z$tm, function(x) {
        max(x) - min(x)
    }, double(1)))
    expect_equal(z$gcContentMean, vapply(z$gcContent, mean, double(1)))
    expect_equal(z$gcContentRange, vapply(z$gcContent, function(x) {
        max(x) - min(x)
    }, double(1)))
    expect_true(all(z$length >= 16))
    expect_true(all(z$length <= 20)) # # end coverage # probe # end runs fwd rev

})

test_that("oligos returns an error when it should", {
    expect_error(oligos(unclass(x)))
    expect_error(oligos(x, maxGapFrequency = -1))
    expect_error(oligos(x, maxGapFrequency = 1.1))
    expect_error(oligos(x, maxGapFrequency = -1))
    expect_error(oligos(x, lengthPrimer = 13))
    expect_error(oligos(x, lengthPrimer = 31))
    expect_error(oligos(x, maxDegeneracyPrimer = 33))
    expect_error(oligos(x, maxDegeneracyPrimer = 0))
    expect_error(oligos(x, gcClampPrimer = 1))
    expect_error(oligos(x, avoidThreeEndRunsPrimer = "true"))
    expect_error(oligos(x, gcRangePrimer = c(-0.1, 1)))
    expect_error(oligos(x, gcRangePrimer = c(0, 1.1)))
    expect_error(oligos(x, tmRangePrimer = c(19, 90)))
    expect_error(oligos(x, tmRangePrimer = c(20, 91)))
    expect_error(oligos(x, concPrimer = 19))
    expect_error(oligos(x, concPrimer = 2001))
    expect_error(oligos(x, designStrategyPrimer = ""))
    expect_error(oligos(x, probe = 1))
    expect_error(oligos(x, lengthProbe = 13))
    expect_error(oligos(x, lengthProbe = 31))
    expect_error(oligos(x, maxDegeneracyProbe = 33))
    expect_error(oligos(x, maxDegeneracyProbe = 0))
    expect_error(oligos(x, maxDegeneracyPrimer = 33))
    expect_error(oligos(x, avoidFiveEndGProbe = 1))
    expect_error(oligos(x, gcRangeProbe = c(-0.1, 1)))
    expect_error(oligos(x, gcRangeProbe = c(0, 1.1)))
    expect_error(oligos(x, tmRangeProbe = c(19, 90)))
    expect_error(oligos(x, tmRangeProbe = c(20, 91)))
    expect_error(oligos(x, concProbe = 19))
    expect_error(oligos(x, concProbe = 2001))
    expect_error(oligos(x, concNa = 0.009))
    expect_error(oligos(x, concNa = 1.1))
})

test_that(".nmers works", {
    seq <- sample(c("A", "C", "G", "T"), 100, replace = TRUE)
    nmer <- .nmers(seq, n = 4)
    expect_equal(ncol(nmer), 4)
    expect_equal(seq[1:4], nmer[1, ])
    expect_equal(seq[2:5], nmer[2, ])
    expect_equal(seq[(length(seq) - 3):length(seq)], nmer[nrow(nmer), ])

    ## Test if it returns a matrix even if it is only one row
    nmer <- .nmers(c("A", "C", "G"), 3)
    expect_true(is.matrix(nmer))
    expect_equal(dim(nmer), c(1, 3))
})

test_that(".countDegeneracy works", {
    seq <- c("A", "C", "R", "N", "Y")
    degen <- .countDegeneracy(seq)
    table <- rprimer:::lookup$degeneracy
    expect_equal(degen, unname(table["R"] * table["N"] * table["Y"]))
    expect_equal(.countDegeneracy(c("A", "C", "G", "T")), 1)
})

test_that(".splitAndPaste works", {
    majority <- .nmers(x$majority, 10)
    iupac <- .nmers(x$iupac, 10)
    fwd <- .splitAndPaste(majority, iupac)
    expect_equal(ncol(fwd), ncol(majority))
    expect_equal(majority[, 1:7], fwd[, 1:7])
    expect_equal(iupac[, 8:10], fwd[, 8:10])
    rev <- .splitAndPaste(iupac, majority, rev = TRUE)
    expect_equal(iupac[, 1:3], rev[, 1:3])
    expect_equal(majority[, 4:10], rev[, 4:10])
    expect_equal(ncol(rev), ncol(majority))
    expect_false(any(grep("[^ACGT-]", fwd[, 1:7])))
    expect_false(any(grep("[^ACGT-]", rev[, 8:10])))
    majority <- majority[1, , drop = FALSE]
    iupac <- iupac[1, , drop = FALSE]
    oneRow <- .splitAndPaste(majority, iupac)
    expect_true(is.matrix(oneRow))
    expect_equal(nrow(oneRow), 1L)
    expect_equal(ncol(oneRow), ncol(majority))
    majority <- .nmers(x$majority, 11)
    iupac <- .nmers(x$iupac, 11)
    fwd <- .splitAndPaste(majority, iupac)
    expect_equal(ncol(fwd), ncol(majority))
    rev <- .splitAndPaste(iupac, majority, rev = TRUE)
    expect_equal(ncol(rev), ncol(majority))
    coverage <- .nmers(x$coverage, 12)
    identity <- .nmers(x$identity, 12)
    identityCoverage <- .splitAndPaste(identity, coverage)
    expect_identical(identityCoverage[, 1:8], identity[, 1:8])
    expect_identical(identityCoverage[, 9:12], coverage[, 9:12])
    coverageIdentity <- .splitAndPaste(coverage, identity, rev = TRUE)
    expect_identical(coverageIdentity[, 1:4], coverage[, 1:4])
    expect_identical(coverageIdentity[, 5:12], identity[, 5:12])
})

test_that(".generateOligos works", {
    expect_equal(
        x$iupac[y$start[4]:y$end[4]], y$iupacSequence[4, ]
    )
    expect_equal(
        y$degeneracy[4],
        .countDegeneracy(x$iupac[y$start[4]:y$end[4]])
    )
    expect_equal(
        y$gapFrequency[4], max(x$gaps[y$start[4]:y$end[4]])
    )
    expect_equal(
        y$coverage[4], mean(x$coverage[y$start[4]:y$end[4]])
    )
    expect_equal(
        y$endCoverageFwd[4],
        min(x$coverage[(y$end[4] - 4):y$end[4]])
    )
    expect_equal(
        y$endCoverageRev[4],
        min(x$coverage[y$start[4]:(y$start[4] + 4)])
    )
    expect_equal(y$roiStart[4], 5000)
    expect_equal(y$roiEnd[4], 6000)
    expect_equal(unique(y$length), ncol(y$iupacSequence))
    expect_true(is.list(y))
})

test_that(".filterOligos works", {
    y <- .filterOligos(y, maxDegeneracy = 32)
    expect_true(all(y$degeneracy <= 32))
    y <- .filterOligos(y, maxGapFrequency = 0)
    expect_true(all(y$gapFrequency == 0))
    expect_true(is.list(y))
})

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

test_that(".makeOligoMatrix works", {
    y <- .generateOligos(x[5000:6000, ])
    y <- .filterOligos(y, maxDegeneracy = 8)
    oligoList <- apply(y$iupacSequence, 1, .expandDegenerates)
    expect_true(all(vapply(oligoList, is.matrix, logical(1))))
    expect_equal(
        ncol(y$iupacSequence), unique(vapply(oligoList, ncol, integer(1)))
    )
    oligoMatrix <- .makeOligoMatrix(oligoList)
    expect_equal(length(unique(rownames(oligoMatrix))), length(oligoList))

    ## Test how it works when the maximum degeneracy is one
    y <- .generateOligos(x[5000:6000, ])
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
    oligoMatrix <- .makeOligoMatrix(oligoList)
    expect_equal(length(unique(rownames(oligoMatrix))), length(oligoList))

    ## Test if there is only one oligo as input
    y <- lapply(y, function(x) {
        if (is.matrix(x)) x[1, , drop = FALSE] else x[1]
    })
    oligoList <- apply(y$iupacSequence, 1, .expandDegenerates)
    if (!is.list(oligoList)) oligoList <- list(t(oligoList))
    oligoMatrix <- .makeOligoMatrix(oligoList)
    expect_equal(length(unique(rownames(oligoMatrix))), length(oligoList))

    ## Test if degeneracy is 32 (max)
    y <- .generateOligos(x[5000:6000, ])
    y <- .filterOligos(y, maxDegeneracy = 32)
    oligoList <- apply(y$iupacSequence, 1, .expandDegenerates)
    expect_true(all(vapply(oligoList, is.matrix, logical(1))))
    expect_equal(
        ncol(y$iupacSequence), unique(vapply(oligoList, ncol, integer(1)))
    )
    oligoMatrix <- .makeOligoMatrix(oligoList)
    expect_equal(length(unique(rownames(oligoMatrix))), length(oligoList))
})

test_that(".reverseComplement works", {
    seq <- t(matrix(c("A", "C", "G")))
    rc <- .reverseComplement(seq)
    expect_true(is.matrix(rc))
    expect_equal(dim(rc), c(1, 3))
    expect_equivalent(rc, c("C", "G", "T"))
})

test_that(".detectGcClamp works", {
    seq <- c("A", "C", "G", "T", "G", "C", "T", "A")
    gc <- as.integer(seq == "C" | seq == "G")
    gc <- t(matrix(gc))
    expect_true(.detectGcClamp(gc))
    seq <- c("A", "C", "G", "T", "T", "C", "T", "A")
    gc <- ifelse(seq == "C" | seq == "G", 1, 0)
    gc <- t(matrix(gc))
    expect_false(.detectGcClamp(gc))
    seq <- c("A", "C", "G", "T", "C", "C", "G", "C")
    gc <- ifelse(seq == "C" | seq == "G", 1, 0)
    gc <- t(matrix(gc))
    expect_false(.detectGcClamp(gc))
    seq <- t(matrix(c("A", "C", "G", "T", "C", "C", "G", "C")))
    rc <- .reverseComplement(seq)
    gcRc <- ifelse(rc == "C" | rc == "G", 1, 0)
    gcRc <- t(matrix(gcRc))
    expect_true(.detectGcClamp(gc, rev = TRUE))
    expect_true(.detectGcClamp(gcRc))
    gc <- matrix(rep(c(1, 0, 1, 1, 0, 0, 0, 0, 0), 8), ncol = 9, byrow = TRUE)
    expect_equal(sum(.detectGcClamp(gc)), 0)
    expect_equal(sum(.detectGcClamp(gc, rev = TRUE)), nrow(gc))
})

test_that(".detectThreeEndRuns works", {
    seq <- t(matrix(c("A", "T", "C", "C", "C")))
    expect_true(.detectThreeEndRuns(seq))
    seq <- t(matrix(c("A", "T", "T", "C", "C")))
    expect_false(.detectThreeEndRuns(seq))
    expect_false(.detectThreeEndRuns(seq, rev = TRUE))
    seq <- t(matrix(c("G", "G", "G", "T", "A")))
    expect_true(.detectThreeEndRuns(seq, rev = TRUE))
    seq <- matrix(rep(c("G", "G", "G", "T", "A"), 10), ncol = 5, byrow = TRUE)
    expect_equal(sum(.detectThreeEndRuns(seq)), 0)
    expect_equal(sum(.detectThreeEndRuns(seq, rev = TRUE)), 10)
})

test_that(".detectRepeats works", {
    mono <- "CTTTTTA"
    di <- "CTCTCTCTCTA"
    expect_true(.detectRepeats(mono))
    expect_true(.detectRepeats(di))
    expect_equal(sum(.detectRepeats(c(di, mono))), 2)
})

test_that(".getAllVariants works", {
    ## Make sure it works if max degeneracy is one
    y <- .filterOligos(y, maxDegeneracy = 1)
    all <- .getAllVariants(y)
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
    y <- .generateOligos(x[100:150, ])
    y <- .filterOligos(y, maxDegeneracy = 16)
    all <- .getAllVariants(y)
    nVariants <- lapply(all, function(x) vapply(x, length, integer(1)))
    nVariants <- do.call("rbind", nVariants)
    expect_true(all(apply(nVariants, 2, function(x) length(unique(x))) == 1))
    expect_true(all(nVariants <= 16))
})

test_that(".getMeanAndRange, and .makeOligoDf work", {
    y <- .generateOligos(x[100:150, ])
    y <- .filterOligos(y, maxDegeneracy = 16)
    all <- .getAllVariants(y)
    meanAndRange <- .getMeanAndRange(all)
    expect_true(is.data.frame(meanAndRange))
    oligoDf <- .makeOligoDf(y)
    expect_true(is.data.frame(oligoDf))
    nrows <- c(nrow(meanAndRange), nrow(oligoDf))
    expect_equal(length(unique(nrows)), 1)

    ## Test with only one entry in each list
    y <- lapply(y, function(x) {
        if (is.matrix(x)) x[1, , drop = FALSE] else x[1]
    })
    oligoDf <- .makeOligoDf(y)
    expect_equal(nrow(oligoDf), 1)
})

test_that(".designOligos and .isWithinRange work", {
    y <- .designOligos(x, maxDegeneracy = 4)
    expect_true(is.data.frame(y))
    expect_true(all(y$degeneracy <= 8))
    expect_error(.designOligos(x[1:50, ], maxDegeneracy = 1))
    gcLogical <- .isWithinRange(y$gcContent, range = c(0.5, 0.55))
    expect_true(is.logical(gcLogical[[1]]))
    gcLocical <- .isWithinRange(y$gcContent[[1]][1], range = c(0.5, 0.55))
    expect_true(is.logical(gcLogical[[1]]))
})

# test_that(".convertToMatrices and .isValid work", {


# })

# test_that(".filterPrimers works", {

# })

# test_that(".filterProbes works", {

# })

# test_that(".beautify works", {

# })
