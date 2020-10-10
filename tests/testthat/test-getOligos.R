test_that(".countDegeneracy works", {
  expect_equal(.countDegeneracy("N"), 4)
  expect_equal(.countDegeneracy("NRC"), 4 * 2 * 1)
  expect_equal(.countDegeneracy("nrc-"), 4 * 2 * 1)
})

test_that(".gcContent works", {
  expect_equal(.gcContent("ACTCTC"), 0.5)
  expect_equal(.gcContent("ACTCTC--"), 0.5)
  expect_equal(.gcContent("actctc--"), 0.5)
})

test_that(".reverseComplement works", {
  expect_equal(.reverseComplement("acg"), "CGT")
  expect_equal(.reverseComplement("ACG"), "CGT")
  expect_equal(.reverseComplement("yrn"), "NYR")
  expect_equal(.reverseComplement("-a"), "T-")
})

test_that(".splitSequence works", {
  expect_equal(.splitSequence("CGGT"), c("C", "G", "G", "T"))
})

test_that(".getNmers works", {
  data("exampleRprimerProperties")
  toTest <- exampleRprimerProperties$Majority[1:50]
  result <- .getNmers(toTest, n = 20)
  expect_true(all(nchar(result) == 20))
})

test_that(".getNmers returns an error when it should", {
  data("exampleRprimerProperties")
  toTest <- exampleRprimerProperties$Majority[1:50]
  expect_error(.getNmers(toTest))
})

test_that(".runningSum works", {
  expect_equal(.runningSum(rep(5, 3), n = 2), c(10, 10))
  expect_equal(.runningSum(rep(1, 10), n = 2), rep(2, 9))
})

test_that(".runningSum returns an error when it should", {
  expect_error(.runningSum(1:10, n = 20))
  expect_error(.runningSum(1:10))
})

test_that(".generateOligos returns an error when it should.", {
  data("exampleRprimerProperties")
  toTest <- exampleRprimerProperties
  expect_error(.generateOligos(toTest, oligoLength = 13))
  expect_error(.generateOligos(toTest, oligoLength = 31))
  expect_error(.generateOligos(toTest, maxGapFrequency = -1))
  expect_error(.generateOligos(toTest, maxGapFrequency = 1.1))
  expect_error(.generateOligos(toTest, maxDegeneracy = 0))
  expect_error(.generateOligos(toTest, maxDegeneracy = 33))
})

test_that(
  ".exclude, .addReverseComplement, .addGcContent, .addTm
  and .filterOligos work", {
  data("exampleRprimerProperties")
  toTest <- .generateOligos(exampleRprimerProperties)
  dinucleotideRepeats <- "(AT){4,}|(TA){4,}|(AC){4,}|(CA){4,}|(AG){4,}|(GA){4,}|(GT){4,}|(TG){4,}|(CG){4,}|(GC){4,}|(TC){4,}|(CT){4,})"
  mononucleotideRepeates <- "([A-Z])\\1\\1\\1\\1"
  expect_true(any(grepl(dinucleotideRepeats, toTest$Majority)))
  expect_true(any(grepl(mononucleotideRepeates, toTest$Majority)))
  toTest <- .exclude(toTest)
  expect_false(any(grepl(dinucleotideRepeats, toTest$Majority)))
  expect_false(any(grepl(mononucleotideRepeates, toTest$Majority)))
  toTest <- .addReverseComplement(toTest)
  expect_equal(
    c("Majority_RC", "IUPAC_RC") %in% colnames(toTest), c(TRUE, TRUE)
  )
  expect_error(.addGcContent(toTest, gcRange = c(0.45, 10)))
  toTest <- .addGcContent(toTest)
  expect_true("GC_majority" %in% colnames(toTest))
  expect_false(any(toTest$GC_majority < 0))
  expect_false(any(toTest$GC_majority > 1))
  expect_error(.addTm(toTest, tmRange =  c(19, 10)))
  expect_error(.addTm(toTest, tmRange =  c(20, 91)))
  toTest <- .addTm(toTest)
  expect_true("Tm_majority" %in% colnames(toTest))
  expect_true(any(grepl("^G", toTest$Majority)))
  expect_true(any(grepl("^G", toTest$Majority_RC)))
  toTest <- .filterOligos(
    toTest, gcClamp = FALSE, avoid5EndG = TRUE, avoid3EndRuns = FALSE)
  expect_false(any(grepl("^G", toTest$Majority)))
  expect_false(any(grepl("^G", toTest$Majority_RC)))
  expect_true(any(grepl("([A-Z])\\1\\1$", toTest$Majority)))
  expect_true(any(grepl("([A-Z])\\1\\1$", toTest$Majority_RC)))
  toTest <- .filterOligos(
    toTest, gcClamp = FALSE, avoid5EndG = TRUE, avoid3EndRuns = TRUE)
  expect_false(any(grepl("([A-Z])\\1\\1$", toTest$Majority)))
  expect_false(any(grepl("([A-Z])\\1\\1$", toTest$Majority_RC)))
})

test_that(".getOligosWithGcClamp works", {
  toTest <- c(
    "CTCTCAAAAA", "CTCTCGGGGG", "CTCTCCCCCC", "CTCTCAGGCA", "CTCTCAGACA",
    "CTCTCAGGGG"
  )
  expect_equal(
    .getOligosWithGcClamp(toTest),
    c(NA, NA, NA, "CTCTCAGGCA", "CTCTCAGACA", NA)
  )
})

test_that(".expandDegenerates works", {
  expect_equal(.expandDegenerates("CAGR"), c("CAGA", "CAGG"))
  expect_equal(.expandDegenerates("CAGN"), c( "CAGA", "CAGC", "CAGG", "CAGT"))
  expect_error(.expandDegenerates("X"))
})

test_that(".nnSplit, .getNnTableValues, .init3End and .init5End work,", {
  expect_equal(.nnSplit("CGGGT"), c("CG", "GG", "GG", "GT"))
  expect_equal(.getNnTableValues("GC", "dH"), -9800)
  expect_equal(.getNnTableValues("GC", "dS"), -24.4)
  expect_equal(.init3End("CGGT"), c("H" = 2300.0, "S" = 4.1))
  expect_equal(.init3End("CGGG"), c("H" = 100.0, "S" = -2.8))
  expect_equal(.init5End("TGGT"), c("H" = 2300.0, "S" = 4.1))
  expect_equal(.init5End("CGGG"), c("H" = 100.0, "S" = -2.8))
})

test_that(".tm works", {
  expect_error(.tm("CTTTGGGGGTTT", concOligo =  0.1e-07))
  expect_error(.tm("CTTTGGGGGTTT", concOligo =  2.1e-06))
  expect_error(.tm("CTTTGGGGGTTT", concNa =  0.0099))
  expect_error(.tm("CTTTGGGGGTTT", concNa =  1.1))
  expect_error(.tm(c("CTTTGGGGGTTT", "CTTTTTTGGGG")))
  tm1 <- .tm("CGGTTTGGC", concOligo = 5e-07)
  tm2 <- .tm("CGGTTTGGC", concOligo = 1e-06)
  expect_true(tm1 < tm2)
  tm3 <- .tm("CGGTTTGGC", concNa = 0.05)
  tm4 <- .tm("CGGTTTGGC", concNa = 0.06)
  expect_true(tm3 < tm4)
})

test_that(".expandOligos works", {
  data("exampleRprimerProperties")
  toTest <- getOligos(
    exampleRprimerProperties[1:100, ], showAllVariants = FALSE
  )
  toTest <- .expandOligos(toTest[1:2, ])
  expect_true(all(c("All", "All_RC", "Tm_all", "GC_all") %in% colnames(toTest)))
  expect_true(is.list(toTest$All[1]))
})
