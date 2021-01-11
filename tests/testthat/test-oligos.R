data("exampleRprimerProfile")
x <- exampleRprimerProfile
oligos <- .generateOligos(x[5000:6000, ])

# test_that(".oligos works", {

# })


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

test_that(".generateOligos works", {
    expect_equal(
        x$iupac[oligos$start[4]:oligos$end[4]], oligos$iupacSequence[4, ]
    )
    expect_equal(
        oligos$degeneracy[4],
        .countDegeneracy(x$iupac[oligos$start[4]:oligos$end[4]])
    )
    expect_equal(
        oligos$gapFrequency[4], max(x$gaps[oligos$start[4]:oligos$end[4]])
    )
    expect_equal(
        oligos$identity[4], mean(x$identity[oligos$start[4]:oligos$end[4]])
    )
    expect_equal(
        oligos$endIdentityFwd[4],
        min(x$identity[(oligos$end[4] - 4):oligos$end[4]])
    )
    expect_equal(
        oligos$endIdentityRev[4],
        min(x$identity[oligos$start[4]:(oligos$start[4] + 4)])
    )
    expect_equal(oligos$roiStart[4], 5000)
    expect_equal(oligos$roiEnd[4], 6000)
    expect_equal(unique(oligos$length), ncol(oligos$iupacSequence))
    expect_true(is.list(oligos))
})

test_that(".filterOligos works", {
    oligos <- .filterOligos(oligos, maxDegeneracy = 32)
    expect_true(all(oligos$degeneracy <= 32))
    oligos <- .filterOligos(oligos, maxGapFrequency = 0)
    expect_true(all(oligos$gapFrequency == 0))
    expect_true(is.list(oligos))
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
    oligos <- .generateOligos(x[5000:6000, ])
    oligos <- .filterOligos(oligos, maxDegeneracy = 8)
    oligoList <- apply(oligos$iupacSequence, 1, .expandDegenerates)
    expect_true(all(vapply(oligoList, is.matrix, logical(1))))
    expect_equal(
        ncol(oligos$iupacSequence), unique(vapply(oligoList, ncol, integer(1)))
    )
    oligoMatrix <- .makeOligoMatrix(oligoList)
    expect_equal(length(unique(rownames(oligoMatrix))), length(oligoList))

    ## Test how it works when the maximum degeneracy is one
    oligos <- .generateOligos(x[5000:6000, ])
    oligos <- .filterOligos(oligos, maxDegeneracy = 1)
    oligoList <- apply(oligos$iupacSequence, 1, .expandDegenerates)
    if (!is.list(oligoList)) {
        oligoList <- t(oligoList)
        oligoList <- lapply(seq_len(nrow(oligoList)), function(i) {
            oligoList[i, , drop = FALSE]
        })
    }
    expect_equal(
        ncol(oligos$iupacSequence), unique(vapply(oligoList, ncol, integer(1)))
    )
    oligoMatrix <- .makeOligoMatrix(oligoList)
    expect_equal(length(unique(rownames(oligoMatrix))), length(oligoList))

    ## Test if there is only one oligo as input
    oligos <- lapply(oligos, function(x) {
        if (is.matrix(x)) x[1, , drop = FALSE] else x[1]
    })
    oligoList <- apply(oligos$iupacSequence, 1, .expandDegenerates)
    if (!is.list(oligoList)) oligoList <- list(t(oligoList))
    oligoMatrix <- .makeOligoMatrix(oligoList)
    expect_equal(length(unique(rownames(oligoMatrix))), length(oligoList))

    ## Test if degeneracy is 32 (max)
    oligos <- .generateOligos(x[5000:6000, ])
    oligos <- .filterOligos(oligos, maxDegeneracy = 32)
    oligoList <- apply(oligos$iupacSequence, 1, .expandDegenerates)
    expect_true(all(vapply(oligoList, is.matrix, logical(1))))
    expect_equal(
        ncol(oligos$iupacSequence), unique(vapply(oligoList, ncol, integer(1)))
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
    gc <- ifelse(seq == "C" | seq == "G", 1, 0)
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
    oligos <- .filterOligos(oligos, maxDegeneracy = 1)
    all <- .getAllVariants(oligos)
    expect_true(is.list(all))
    expect_equal(
        all$sequence[[1]], paste(oligos$iupacSequence[1, ], collapse = "")
    )
    expect_equal(
        all$sequence[[length(all$sequence)]],
        paste(oligos$iupacSequence[nrow(oligos$iupacSequence), ], collapse = "")
    )
    expect_equal(unique(vapply(all$sequence, length, integer(1))), 1)

    ## If degeneracy is high
    oligos <- .generateOligos(x[100:150, ])
    oligos <- .filterOligos(oligos, maxDegeneracy = 16)
    all <- .getAllVariants(oligos)
    nVariants <- lapply(all, function(x) vapply(x, length, integer(1)))
    nVariants <- do.call("rbind", nVariants)
    expect_true(all(apply(nVariants, 2, function(x) length(unique(x))) == 1))
    expect_true(all(nVariants <= 16))
})

test_that(".getMeanAndRange, and .makeOligoDf work", {
    oligos <- .generateOligos(x[100:150, ])
    oligos <- .filterOligos(oligos, maxDegeneracy = 16)
    all <- .getAllVariants(oligos)
    meanAndRange <- .getMeanAndRange(all)
    expect_true(is.data.frame(meanAndRange))
    #   expect_equal(meanAndRange$tmMean[1], mean(all$tm[[1]]))
    #   expect_equal(meanAndRange$tmRange[1], max(all$tm[[1]]) - min(all$tm[[1]]))
    oligoDf <- .makeOligoDf(oligos)
    expect_true(is.data.frame(oligoDf))
    nrows <- c(nrow(meanAndRange), nrow(oligoDf))
    expect_equal(length(unique(nrows)), 1)

    ## Test with only one entry in each list
    oligos <- lapply(oligos, function(x) {
        if (is.matrix(x)) x[1, , drop = FALSE] else x[1]
    })
    oligoDf <- .makeOligoDf(oligos)
    expect_equal(nrow(oligoDf), 1)
})

test_that(".designOligos and .isWithinRange work", {
    oligos <- .designOligos(x, maxDegeneracy = 4)
    expect_true(is.data.frame(oligos))
    expect_true(all(oligos$degeneracy <= 8))
    expect_error(.designOligos(x[1:50, ], maxDegeneracy = 1))
    gcLogical <- .isWithinRange(oligos$gcContent, range = c(0.5, 0.55))
    expect_true(is.logical(gcLogical[[1]]))
    gcLocical <- .isWithinRange(oligos$gcContent[[1]][1], range = c(0.5, 0.55))
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
