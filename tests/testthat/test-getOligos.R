data("exampleRprimerProfile")
testdata <- exampleRprimerProfile[1:500, ]

test_that("getOligos works", {
    oligos <- getOligos(testdata,
        lengthPrimer = 18:20,
        maxDegeneracyPrimer = 8,
        gcRangePrimer = c(0.50, 0.55),
        tmRangePrimer = c(55, 60),
        probe = FALSE
    )
    expect_s4_class(oligos, "RprimerOligo")
    expect_false(any(oligos$type == "probe"))
    expect_false(any(oligos$length < 18))
    expect_false(any(oligos$length > 20))
    expect_false(any(oligos$degeneracy > 8))
    expect_false(any(oligos$gcMajority < 0.50))
    expect_false(any(oligos$gcMajority > 0.55))
    expect_false(any(oligos$tmMajority < 55))
    expect_false(any(oligos$tmMajority > 60))
    oligos <- getOligos(testdata, probe = TRUE)
    expect_true(any(oligos$type == "probe"))
    oligos <- oligos[oligos$type == "probe", ]
    expect_false(any(oligos$length < 18))
    expect_false(any(oligos$length > 22))
    expect_false(any(oligos$degeneracy > 4))
    expect_false(any(oligos$gcMajority < 0.45))
    expect_false(any(oligos$gcMajority > 0.55))
    expect_false(any(oligos$tmMajority < 55))
    expect_false(any(oligos$tmMajority > 70))
})

test_that("getOligos returns an error when it should", {
    expect_error(getOligos(unclass(testdata)))
    expect_error(getOligos(testdata, showAllVariants = "TRUE"))
    expect_error(getOligos(testdata[1:20, ]))
    expect_error(getOligos(testdata[1:100, ], maxDegeneracyPrimer = 1))
    expect_error(getOligos(testdata[1:100, ], maxDegeneracyProbe = 1))
})

test_that(".getNmers works", {
    nmers <- .getNmers(testdata$majority, n = 5)
    expect_equal(nmers[1], paste(testdata$majority[1:5], collapse = ""))
    expect_equal(nmers[2], paste(testdata$majority[2:6], collapse = ""))
    expect_true(all(nchar(nmers) == 5))
    expect_equal(
        nmers[length(nmers)], paste(testdata$majority[
            (length(testdata$majority) - 4):length(testdata$majority)
        ], collapse = "")
    )
})

test_that(".getNmers returns an error when it should", {
    expect_error(nmers <- .getNmers(testdata$majority[1:20], n = 22))
})

test_that(".splitSequence works", {
    expect_equal(.splitSequence("ACCG"), c("A", "C", "C", "G"))
})

test_that(".runningSum works", {
    expect_equal(.runningSum(rep(5, 3), n = 2), c(10, 10))
    expect_equal(.runningSum(rep(1, 10), n = 2), rep(2, 9))
})

test_that(".runningSum returns an error when it should", {
    expect_error(.runningSum(1:10, n = 20))
    expect_error(.runningSum(1:10))
})

test_that(".countEndIdentity works", {
    identity <- testdata$identity[100:119]
    endIdentity <- .countEndIdentity(identity, n = 19)
    expect_equivalent(
        endIdentity[1, 1],
        min(identity[(length(identity) - 5):(length(identity) - 1)])
    )
    expect_equivalent(
        endIdentity[2, 1],
        min(identity[(length(identity) - 4):length(identity)])
    )
    expect_equivalent(endIdentity[1, 2], min(identity[1:5]))
    expect_equivalent(endIdentity[2, 2], min(identity[2:6]))
})

test_that(".countDegeneracy works", {
    expect_equal(.countDegeneracy("N"), 4)
    expect_equal(.countDegeneracy("NRC"), 4 * 2 * 1)
    expect_equal(.countDegeneracy("nrc-"), 4 * 2 * 1)
    expect_true(is.na(.countDegeneracy("X")))
})

test_that(".generateOligos works", {
    oligos <- .generateOligos(
        testdata,
        oligoLength = 20,
        maxGapFrequency = 0.05, maxDegeneracy = 8
    )
    expect_true(all(nchar(oligos$majority) == 20))
    expect_true(all(oligos$degeneracy <= 8))
    invalidPositions <- which(testdata$gaps > 0.05)
    oligoPositions <- purrr::map2(oligos$start, oligos$end, function(x, y) {
        x:y
    })
    oligoPositions <- do.call("rbind", oligoPositions)
    expect_equal(length(intersect(invalidPositions, oligoPositions)), 0)
})

test_that(".generateOligos returns an error when it should", {
    expect_error(.generateOligos(testdata, oligoLength = 13))
    expect_error(.generateOligos(testdata, oligoLength = 31))
    expect_error(.generateOligos(testdata, maxGapFrequency = -0.1))
    expect_error(.generateOligos(testdata, maxGapFrequency = 1.1))
    expect_error(.generateOligos(testdata, maxDegeneracy = 0))
    expect_error(.generateOligos(testdata, maxDegeneracy = 33))
})

test_that(".exclude works", {
    oligos <- data.frame( # Two invalid oligos and two valid
        majority = c("ATATATAT", "ATCCCCCT", "ATATATCC", "ATCCCCT"),
        mockvector = rep(NA, 4)
    )
    dinucleotideRepeats <- "(AT){4,}|(TA){4,}|(AC){4,}|(CA){4,}|(AG){4,}|(GA){4,}|(GT){4,}|(TG){4,}|(CG){4,}|(GC){4,}|(TC){4,}|(CT){4,})"
    expect_true(any(grepl(dinucleotideRepeats, oligos$majority)))
    mononucleotideRepeates <- "([A-Z])\\1\\1\\1\\1"
    expect_true(any(grepl(mononucleotideRepeates, oligos$majority)))
    oligos <- .exclude(oligos)
    expect_false(any(grepl(dinucleotideRepeats, oligos$majority)))
    expect_false(any(grepl(mononucleotideRepeates, oligos$majority)))
    expect_equal(nrow(oligos), 2)
})

test_that(".reverseComplement works", {
    expect_equal(.reverseComplement("acg"), "CGT")
    expect_equal(.reverseComplement("ACG"), "CGT")
    expect_equal(.reverseComplement("yrn"), "NYR")
    expect_equal(.reverseComplement("-a"), "T-")
})

test_that(".addReverseComplement works", {
    oligos <- .generateOligos(testdata[1:50, ])
    # Should work even if there is 0 oligos
    oligosEmpty <- oligos[-(1:nrow(oligos)), ]
    oligos <- .addReverseComplement(oligos)
    expect_true(all(c("majorityRc", "iupacRc") %in% colnames(oligos)))
    oligosEmpty <- .addReverseComplement(oligosEmpty)
    expect_true(all(c("majorityRc", "iupacRc") %in% colnames(oligosEmpty)))
})

test_that(".gcContent works", {
    expect_equal(.gcContent("ACTCTC"), 0.5)
    expect_equal(.gcContent("ACTCTC--"), 0.5)
    expect_equal(.gcContent("actctc--"), 0.5)
})

test_that(".addGcContent work", {
    oligos <- .generateOligos(testdata[1:50, ])
    # Should work even if there is 0 oligos
    oligosEmpty <- oligos[-(1:nrow(oligos)), ]
    oligos <- .addGcContent(oligos)
    expect_true("gcMajority" %in% colnames(oligos))
    oligosEmpty <- .addGcContent(oligosEmpty)
    expect_true("gcMajority" %in% colnames(oligosEmpty))
})

test_that(".addGcContent returns an error when it should", {
    oligos <- .generateOligos(testdata[1:50, ])
    expect_error(.addGcContent(oligos, gcRange = c(-0.1, 1)))
    expect_error(.addGcContent(oligos, gcRange = c(0, 1.1)))
})

test_that(".getOligosWithGcClamp works", {
    oligos <- c(
        "CTCTCAAAAA", "CTCTCGGGGG", "CTCTCCCCCC", "CTCTCAGGCA", "CTCTCAGACA",
        "CTCTCAGGGG"
    )
    expect_equal(
        .getOligosWithGcClamp(oligos),
        c(NA, NA, NA, "CTCTCAGGCA", "CTCTCAGACA", NA)
    )
})

test_that(".filterOligos work", {
    oligos <- data.frame(
        majority = c(
            "GCACGCGCGGT", # 5 end G
            "CGTCCTTAATA", # no GC clamp
            "CTTTCCCAAAA" # 3 end runs
        ),
        majorityRc = c(
            .reverseComplement("GCACGCGCGGT"),
            .reverseComplement("CGTCCTTAATA"),
            .reverseComplement("CTTTCCCAAAA")
        )
    )
    expect_true(any(grepl("^G", oligos$majority)))
    expect_true(any(grepl("([A-Z])\\1\\1$", oligos$majority)))
    expect_false(is.na(oligos$majority[2]))
    oligos <- .filterOligos(
        oligos,
        gcClamp = TRUE, avoid5EndG = TRUE, avoid3EndRuns = TRUE,
        minEndIdentity = NULL
    )
})

test_that(".filterOligos returns an error when it should.", {
    oligos <- .generateOligos(testdata[1:50, ])
    expect_error(.filterOligos(oligos, gcClamp = "TRUE"))
    expect_error(.filterOligos(oligos, avoid5EndG = "TRUE"))
    expect_error(.filterOligos(oligos, avoid3EndRuns = "TRUE"))
    expect_error(.filterOligos(oligos, minEndIdentity = 1.1))
    expect_error(.filterOligos(oligos, minEndIdentity = -0.1))
})

test_that(".nnSplit, .getNnTableValues, .init3End and .init5End work,", {
    expect_equal(.nnSplit("CGGGT"), c("CG", "GG", "GG", "GT"))
    expect_equal(.getStackValue("GC", "dH"), -9800)
    expect_equal(.getStackValue("GC", "dS"), -24.4)
})

test_that(".tm works", {
    expect_error(.tm("CTTTGGGGGTTT", concOligo = 19))
    expect_error(.tm("CTTTGGGGGTTT", concOligo = 2001))
    expect_error(.tm("CTTTGGGGGTTT", concNa = 0.0099))
    expect_error(.tm("CTTTGGGGGTTT", concNa = 1.1))
    expect_error(.tm(c("CTTTGGGGGTTT", "CTTTTTTGGGG")))
    tm1 <- .tm("CGGTTTGGC", concOligo = 500)
    tm2 <- .tm("CGGTTTGGC", concOligo = 1000)
    expect_true(tm1 < tm2)
    tm3 <- .tm("CGGTTTGGC", concNa = 0.05)
    tm4 <- .tm("CGGTTTGGC", concNa = 0.06)
    expect_true(tm3 < tm4)
})

test_that(".addTm works", {
    oligos <- .generateOligos(testdata)
    oligos <- .addTm(oligos)
    expect_true("tmMajority" %in% names(oligos))
    oligos <- .generateOligos(testdata)
    oligos <- .addTm(oligos[-(1:nrow(oligos)), ])
    expect_true("tmMajority" %in% names(oligos))
})

test_that(".addTm returns an error when it should", {
    oligos <- .generateOligos(testdata)
    expect_error(.addTm(oligos, tmRange = c(19, 50)))
    expect_error(.addTm(oligos, tmRange = c(50, 91)))
})

test_that(".expandDegenerates works", {
    expect_equal(.expandDegenerates("CAGR"), c("CAGA", "CAGG"))
    expect_equal(.expandDegenerates("CAGN"), c("CAGA", "CAGC", "CAGG", "CAGT"))
    expect_error(.expandDegenerates("X"))
})
