## Load example data
data("exampleRprimerAlignment")
data("exampleRprimerOligo")
data("exampleRprimerAssay")
data("exampleRprimerMatchOligo")
data("exampleRprimerMatchAssay")

Biostrings::rowmask(exampleRprimerAlignment, invert = TRUE) <- 1

# checkMatch ===================================================================

test_that("checkMatch returns an error when it should", {
    expect_error(
        checkMatch(exampleRprimerOligo[1:2, ], unclass(exampleRprimerAlignment))
    )
    expect_error(
        checkMatch(exampleRprimerAssay[1:2, ], unclass(exampleRprimerAlignment))
    )
})

test_that("checkMatch works", {
    oligoMatch <- checkMatch(exampleRprimerOligo[1, ], exampleRprimerAlignment)
    assayMatch <- checkMatch(exampleRprimerAssay[1, ], exampleRprimerAlignment)
    expect_s4_class(oligoMatch, "RprimerMatchOligo")
    expect_s4_class(assayMatch, "RprimerMatchAssay")
})

# .maskRange ===================================================================

# .getMatchIndex ===============================================================

test_that(".getMatchIndex works", {
    target <- c(
        "PM" = "AAAAAA",
        "MM1" = "AAAAAT",
        "MM2" = "AAAATT",
        "MM3" = "AAATTT",
        "MM4" = "AATTTT"
    )
    target <- Biostrings::DNAStringSet(target)
    x <- Biostrings::DNAStringSet("AAAAAA") ## Oligo
    result <- .getMatchIndex(x, target)
    expectedResult <- list(1, 2, 3, 4, 5)
    expect_equal(result, expectedResult)
    xDegenerate <- Biostrings::DNAStringSet(c("AAAAAA", "AAAAAT"))
    result <- .getMatchIndex(xDegenerate, target)
    expectedResult <- list(c(1, 2), 3, 4, 5, integer(0L))
    expect_equal(result, expectedResult)
})

# .getSequenceNames ============================================================

test_that(".getSequenceNames works", {
    target <- list("A BLABLA" = "CCC", "B BLABLA" = "TTT")
    result <- .getSequenceNames(1, target)
    expect_equal(result, list("A"))
})

# .getMatchProportion ==========================================================

test_that(".getMatchProportion works", {
    x <- Biostrings::DNAStringSet("AAAAA")
    target <- c(
        "AAAAA",
        "AAAAT",
        "AAATT",
        "AATTT",
        "ATTTT",
        "TTTTT"
    )
    target <- Biostrings::DNAStringSet(target)
    result <- .getMatchProportion(x, target)
    expect_true(is.data.frame(result))
})

# .checkMatchOligo =============================================================

test_that(".checkMatchOligo works", {
    x <- as.data.frame(exampleRprimerOligo[1, ])

    ## Check that target and off target identification works
    onTarget <- exampleRprimerAlignment
    Biostrings::colmask(onTarget, invert = TRUE) <- x$start:x$end
    on <- .checkMatchOligo(x, onTarget)
    offTarget <- exampleRprimerAlignment
    Biostrings::colmask(offTarget) <- x$start:x$end
    off <- .checkMatchOligo(x, offTarget)
    both <- .checkMatchOligo(x, exampleRprimerAlignment)
    expect_equal(on[!grepl("Off", names(on))], both[!grepl("Off", names(both))])
    expect_equal(off[grepl("Off", names(off))], both[grepl("Off", names(both))])
})

# .convertAssay and .checkMatchAssay ===========================================

test_that(".convertAssay and checkMatchAssay work", {
    x <- as.data.frame(exampleRprimerAssay[1:2, ])
    x <- x[!grepl("Pr$", names(x))]
    result <- .convertAssay(x)
    expect_true(is.list(result))
    result <- .checkMatchAssay(x, exampleRprimerAlignment)
    expect_true(is.data.frame(result))
})
