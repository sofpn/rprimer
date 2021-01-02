## Data
data("exampleRprimerProfile")
x <- exampleRprimerProfile
oligos <- .generateOligos(x[5000:6000, ])

test_that(".getNmers works", {
  seq <- sample(c("A", "C", "G", "T"), 100, replace = TRUE)
  nmer <- .getNmers(seq, n = 4)
  expect_equal(ncol(nmer), 4)
  expect_equal(seq[1:4], nmer[1, ])
  expect_equal(seq[2:5], nmer[2, ])
  expect_equal(seq[(length(seq) - 3):length(seq)], nmer[nrow(nmer), ])

  ## Test if it returns a matrix even if it is only one row
  nmer <- .getNmers(c("A", "C", "G"), 3)
  expect_true(is.matrix(nmer))
  expect_equal(dim(nmer), c(1, 3))
})

test_that(".countDegeneracy works", {
  seq <- c("A", "C", "R", "N", "Y")
  degen <- .countDegeneracy(seq)
  table <- rprimer:::lookup$degeneracy
  expect_equal(degen, unname(table["R"]*table["N"]*table["Y"]))
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
  expect_equal(oligos$alignmentStart[4], 400)
  expect_equal(oligos$alignmentEnd[4], 600)
  expect_equal(unique(oligos$length), ncol(oligos$iupacSequence))
  expect_true(is.list(oligos))
})

test_that(".filterOligos works", {
  oligos <- .filterOligos(oligos, maxDegeneracy = 32)
  expect_true(all(oligos$degeneracy <= 32))
  oligos <- .filterOligos(oligos, maxGapFrequency = 0)
  expect_true(all(oligos$degeneracy == 0))
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
  if (!is.list(oligoList)) oligoList <- list(t(oligoList))
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
  seq <- c("A", "C", "G")
  rc <- .reverseComplement(seq)
  expect_true(is.matrix(rc))
  expect_equal(dim(rc), c(1, 3))
  expect_equivalent(rc, c("C", "G", "T"))
})

test_that(".detectGcClamp works", {
  seq <- c("A", "C", "G", "T", "G", "C", "T", "A")
  gc <- ifelse(seq == "C" | seq == "G", 1, 0)
  expect_true(.detectGcClamp(gc))
  seq <- c("A", "C", "G", "T", "T", "C", "T", "A")
  gc <- ifelse(seq == "C" | seq == "G", 1, 0)
  expect_false(.detectGcClamp(gc))
  seq <- c("A", "C", "G", "T", "C", "C", "G", "C")
  gc <- ifelse(seq == "C" | seq == "G", 1, 0)
  expect_false(.detectGcClamp(gc))
  seq <- c("A", "C", "G", "T", "C", "C", "G", "C")
  rc <- as.vector(.reverseComplement(seq))
  gcRc <- ifelse(rc == "C" | rc == "G", 1, 0)
  expect_true(.detectGcClamp(gc, fwd = FALSE))
  expect_true(.detectGcClamp(gcRc))
  gc <- matrix(rep(c(1, 0, 1, 1, 0, 0, 0, 0, 0), 8), ncol = 9, byrow = TRUE)
  expect_equal(sum(.detectGcClamp(gc)), 0)
  expect_equal(sum(.detectGcClamp(gc, fwd = FALSE)), nrow(gc))
})

test_that(".detectThreeEndRuns works", {
  seq <- c("A", "T", "C", "C", "C")
  expect_true(.detectThreeEndRuns(seq))
  seq <- c("A", "T", "T", "C", "C")
  expect_false(.detectThreeEndRuns(seq))
  expect_false(.detectThreeEndRuns(seq, fwd = FALSE))
  seq <- c("G", "G", "G", "T", "A")
  expect_true(.detectThreeEndRuns(seq, fwd = FALSE))
  seq <- matrix(rep(c("G", "G", "G", "T", "A"), 10), ncol = 5, byrow = TRUE)
  expect_equal(sum(.detectThreeEndRuns(seq)), 0)
  expect_equal(sum(.detectThreeEndRuns(seq, fwd = FALSE)), 10)
})

test_that(".getAllVariants works", {

  ## Make sure it works if max degeneracy is one
  x <- .filterOligos(oligos, maxDegeneracy = 1, maxGapFrequency = 0.1)
  all <- .getAllVariants(x)

})

