#' Design oligos
#'
#' \code{oligos()} designs oligos (primers and probes)
#' from an \code{RprimerProfile} object.
#'
#' @param x An \code{RprimerProfile} object.
#'
#' @param maxGapFrequency
#' Maximum allowed gap frequency (for both primers and probes).
#' A number [0, 1], defaults to 0.05.
#'
#' @param lengthPrimer
#' Primer length. A numeric vector [15, 30],
#' defaults to \code{18:22}.
#'
#' @param maxDegeneracyPrimer
#' Maximum number of variants of each primer. A number [1, 64], defaults to 4.
#'
#' @param gcClampPrimer
#' If primers with no GC-clamp
#' should be avoided.
#' A GC-clamp
#' is identified as two to three G or
#' C:s within the last five bases (3'-end) of the oligo.
#' \code{TRUE} or \code{FALSE}, defaults to \code{TRUE}.
#'
#' @param avoidThreeEndRunsPrimer
#' If primers with more than two runs
#' of the same nucleotide at the terminal 3'-end should be avoided.
#' \code{TRUE} or \code{FALSE}, defaults to \code{TRUE}.
#'
#' @param minThreeEndCoveragePrimer
#' Minimum allowed coverage at the 3' end (the last five bases).
#' A number [0, 1]. If set to 1, all bases of the 3' end must cover all sequence
#' variants in the target alignment. Defaults to 0.98.
#'
#' @param gcRangePrimer
#' GC-content range for primers.
#' A numeric vector [0, 1], defaults to \code{c(0.40, 0.65)}.
#'
#' @param tmRangePrimer
#' Tm range for primers.
#' A numeric vector [30, 90], defaults to \code{c(55, 65)}.
#'
#' @param concPrimer
#' Primer concentration in nM, for Tm calculation. A number
#' [20, 2000], defaults to 500.
#'
#' @param designStrategyPrimer
#' \code{"ambiguous"} or \code{"mixed"}. Defaults to \code{"ambiguous"}.
#'
#' @param probe
#' If probes should be designed. \code{TRUE} or \code{FALSE},
#' defaults to \code{TRUE}.
#'
#' @param lengthProbe
#' Probe length. A numeric vector [15, 40],
#' defaults to \code{18:22}.
#'
#' @param maxDegeneracyProbe
#' Maximum number of variants of each probe. A number [1, 64], defaults to 4.
#'
#' @param avoidFiveEndGProbe
#' If probes with G
#' at the 5'-end should be avoided. \code{TRUE} or \code{FALSE},
#' defaults to \code{TRUE}.
#'
#' @param gcRangeProbe
#' GC-content range for probes (proportion). A numeric vector [0, 1],
#' defaults to \code{c(0.40, 0.65)}.
#'
#' @param tmRangeProbe
#' Tm range for probes.
#' A numeric vector [30, 90], defaults to \code{c(55, 70)}.
#'
#' @param concProbe
#' Primer concentration in nM, for Tm calculation. A numeric vector
#' [20, 2000], defaults to 250.
#'
#' @param concNa
#' The sodium ion concentration in the PCR reaction in M, for Tm calculation.
#' A numeric vector [0.01, 1], defaults to 0.05 (50 mM).
#'
#' @section Output:
#'
#' The output contains the following information:
#'
#' \describe{
#'   \item{type}{Whether the oligo is a primer or probe.}
#'   \item{fwd}{\code{TRUE} if the oligo is valid in forward direction,
#'     \code{FALSE} otherwise.}
#'   \item{rev}{\code{TRUE} if the oligo is valid in reverse direction,
#'     \code{FALSE} otherwise.}
#'   \item{start}{Start position of the oligo.}
#'   \item{end}{End positon of the oligo.}
#'   \item{length}{Oligo length.}
#'   \item{iupacSequence}{Oligo sequence, with amiguous bases (if any).}
#'   \item{iupaSequenceRc}{The reverse complement of the iupacSequence.}
#'   \item{identity}{For ambiguous oligos: Average identity of the oligo.
#'     For mixed oligos: Average identity of the 5' (consensus) part of the
#'     oligo. The value can range from 0 to 1.}
#'   \item{coverage}{For ambiguous oligos: Average coverage of the oligo.
#'     For mixed oligos: Average coverage of the 3' (degenerate) part of the
#'     oligo. The value can range from 0 to 1.}
#'   \item{degeneracy}{Number of sequence variants of the oligo.}
#'   \item{gcContentMean}{Mean GC-content of all sequence variants of the oligo.
#'   }
#'   \item{gcContentRange}{Range in GC-content of all sequence variants of
#'     the oligo.}
#'   \item{tmMean}{Mean tm of all sequence variants of the oligo.}
#'   \item{tmRange}{Range in tm of all sequence variants of the oligo.}
#'   \item{sequence}{All sequence variants of the oligo.}
#'   \item{sequenceRc}{Reverse complements of all sequence variants.}
#'   \item{gcContent}{GC-content of all sequence variants.}
#'   \item{tm}{Tm of all sequence variants.}
#'   \item{method}{Design method used to generate the oligo: "ambiguous",
#'   "mixedFwd" or "mixedRev".}
#'   \item{roiStart}{First position of the input \code{RprimerProfile} object
#'     (roi = region of interest).}
#'   \item{roiEnd}{Last position of the input \code{RprimerProfile} object.}
#' }
#'
#' @section Design strategy for primers:
#'
#' Primers can be designed in either one of two ways: "ambiguous" or "mixed".
#'
#' \itemize{
#' \item{The ambiguous strategy (default) generates primers from the IUPAC
#' consensus sequence, meaning that degenerate bases can occur at
#' any position in the oligo. Probes are always designed using the ambiguous
#' strategy}
#' \item{The mixed strategy is partially based on the Consensus-Degenerate
#' Hybrid
#' Oligonucleotide Primer (CODEHOP) principle (Rose et al. 1998).
#' It is recommended for highly
#' variable targets. Here, primers are generated from both the
#' majority and the IUPAC
#' consensus sequence, and consist of a shorter degenerate part at the
#' 3' end (~1/3 of the
#' primer, targeting a conserved region), and a longer consensus part
#' at the 5' end (~2/3 of the primer)}
#' }
#'
#' @section Validity checks:
#'
#' For an oligo to be considered as valid, all sequence variants must fulfill
#' all the specified design constraints.
#'
#' Oligos with sequence variants containing
#' more than four consecutive runs
#' of the same
#' nucleotide (e.g. "AAAAA") and/or more than three consecutive runs
#' of the same di-nucleotide (e.g. "TATATATA") are excluded.
#'
#' @section Tm-calculation:
#'
#' Melting temperatures are calculated using SantaLucia's nearest-neighbor
#' method, with the following assumptions:
#'
#' \itemize{
#'   \item Oligos are not expected to be self-complementary (no symmetry
#'   correction is done).
#'   \item The oligo concentration is assumed to be much higher
#'   than the target concentration.
#' }
#'
#' See references for table values and equations. Table values can also be
#' found by calling \code{rprimer:::lookup$nn}.
#'
#' @return
#' An \code{RprimerOligo} object. An error message will return
#' if no oligos are found.
#'
#' @references
#' Rose, Timothy M, Emily R Schultz, Jorja G Henikoff, Shmuel Pietrokovski,
#' Claire M McCallum, and Steven Henikoff. 1998. "Consensus-Degenerate
#' Hybrid Oligonucleotide Primers for
#' Amplification of Distantly Related Sequences." Nucleic Acids Research 26 (7):
#' 1628-35.
#'
#' SantaLucia Jr, J., & Hicks, D. (2004).
#' The thermodynamics of DNA structural motifs.
#' Annu. Rev. Biophys. Biomol. Struct., 33, 415-440.
#'
#' @export
#'
#' @examples
#' data("exampleRprimerProfile")
#' target <- exampleRprimerProfile
#'
#' ## Design primers and probes with default values
#' oligos(target)
#'
#' ## Design primers and probes only within a specific region of interest
#' roi <- target[target$position >= 5000 & target$position <= 6000, ]
#' oligos(roi)
#'
#' ## Design primers only
#' oligos(roi, probe = FALSE)
#'
#' ## Allow higher degeneracy
#' oligos(roi,
#'     maxDegeneracyPrimer = 32,
#'     probe = FALSE
#' )
oligos <- function(x,
                   maxGapFrequency = 0.05,
                   lengthPrimer = 18:22,
                   maxDegeneracyPrimer = 4,
                   gcClampPrimer = TRUE,
                   avoidThreeEndRunsPrimer = TRUE,
                   minThreeEndCoveragePrimer = 0.98,
                   gcRangePrimer = c(0.40, 0.65),
                   tmRangePrimer = c(50, 65),
                   concPrimer = 500,
                   designStrategyPrimer = "ambiguous",
                   probe = TRUE,
                   lengthProbe = 18:22,
                   maxDegeneracyProbe = 4,
                   avoidFiveEndGProbe = TRUE,
                   gcRangeProbe = c(0.40, 0.65),
                   tmRangeProbe = c(50, 70),
                   concProbe = 250,
                   concNa = 0.05) {
    if (!methods::is(x, "RprimerProfile")) {
        stop("'x' must be an RprimerProfile object.", call. = FALSE)
    }
    if (!(maxGapFrequency >= 0 && maxGapFrequency <= 1)) {
        stop("'maxGapFrequency' must be from 0 to 1.", call. = FALSE)
    }
    if (!(min(lengthPrimer) >= 15 && max(lengthPrimer) <= 40)) {
        stop("'lengthPrimer' must be from 15 to 40.", call. = FALSE)
    }
    if (!(maxDegeneracyPrimer >= 1 && maxDegeneracyPrimer <= 64)) {
        stop("'maxDegeneracyPrimer' must be from 1 to 64.", call. = FALSE)
    }
    if (!is.logical(gcClampPrimer)) {
        stop("'gcClampPrimer' must be TRUE or FALSE", call. = FALSE)
    }
    if (!is.logical(avoidThreeEndRunsPrimer)) {
        stop("'avoidThreeEndRunsPrimer' must be TRUE or FALSE", call. = FALSE)
    }
    if (!(minThreeEndCoveragePrimer >= 0 && minThreeEndCoveragePrimer <= 1)) {
        stop("'minThreeEndCoveragePrimer' must be from 0 to 1", call. = FALSE)
    }
    if (!(min(gcRangePrimer) >= 0 && max(gcRangePrimer) <= 1)) {
        stop(
            "'gcRangePrimer' must be from 0 to 1, e.g. c(0.45, 0.65).",
            call. = FALSE
        )
    }
    if (!(min(tmRangePrimer) >= 20 && max(tmRangePrimer) <= 90)) {
        stop(
            "'tmRangePrimer' must be from 20 to 90, e.g. c(55, 60).",
            call. = FALSE
        )
    }
    if (!(concPrimer >= 20 && concPrimer <= 2000)) {
        stop("'concPrimer' must be from 20 to 2000.", call. = FALSE)
    }
    if (!(
        designStrategyPrimer == "ambiguous" || designStrategyPrimer == "mixed")
    ) {
        stop(
            "'designStrategyPrimer' must be either 'ambiguous' or 'mixed'.",
            call. = FALSE
        )
    }
    if (!is.logical(probe)) {
        stop("'probe' must be TRUE or FALSE", call. = FALSE)
    }
    if (!(min(lengthProbe) >= 15 && max(lengthProbe) <= 40)) {
        stop("'lengthProbe' must be from 15 to 40.", call. = FALSE)
    }
    if (!(maxDegeneracyProbe >= 1 && maxDegeneracyProbe <= 64)) {
        stop("'maxDegeneracyProbe' must be from 1 to 64.", call. = FALSE)
    }
    if (!is.logical(avoidFiveEndGProbe)) {
        stop("'avoidFiveEndGProbe' must be TRUE or FALSE", call. = FALSE)
    }
    if (!(min(gcRangeProbe) >= 0 && max(gcRangeProbe) <= 1)) {
        stop(
            "'gcRangeProbe' must be from 0 to 1, e.g. c(0.45, 0.65).",
            call. = FALSE
        )
    }
    if (!(min(tmRangeProbe) >= 20 && max(tmRangeProbe) <= 90)) {
        stop(
            "'tmRangeProbe' must be from 20 to 90, e.g. c(55, 60).",
            call. = FALSE
        )
    }
    if (!(concProbe >= 20 && concProbe <= 2000)) {
        stop("'concProbe' must be from 20 to 2000.", call. = FALSE)
    }
    if (!(concNa >= 0.01 && concNa <= 1)) {
        stop("'concNa' must be from 0.01 to 1.", call. = FALSE)
    }
    lengthOligo <- lengthPrimer
    if (probe) {
        lengthOligo <- unique(c(lengthOligo, lengthProbe))
    }
    lengthOligo <- lengthOligo[order(lengthOligo)]
    maxDegeneracy <- maxDegeneracyPrimer
    if (probe) {
        maxDegeneracy <- max(c(maxDegeneracyPrimer, maxDegeneracyProbe))
    }
    oligos <- .designOligos(
        x,
        lengthOligo,
        maxGapFrequency,
        maxDegeneracy,
        concPrimer,
        designStrategyPrimer,
        probe,
        concProbe,
        concNa
    )
    primers <- .filterPrimers(oligos,
        lengthPrimer,
        maxDegeneracyPrimer,
        gcClampPrimer,
        avoidThreeEndRunsPrimer,
        minThreeEndCoveragePrimer,
        gcRangePrimer,
        tmRangePrimer,
        designStrategyPrimer,
        rowThreshold = 1,
        colThreshold = 1
    )
    if (nrow(primers) == 0L) {
        stop("No primers were found.", call. = FALSE)
    }
    if (probe) {
        probes <- .filterProbes(oligos,
            lengthProbe,
            maxDegeneracyProbe,
            avoidFiveEndGProbe,
            gcRangeProbe,
            tmRangeProbe,
            rowThreshold = 1,
            colThreshold = 1
        )
        if (nrow(probes) == 0L) {
            stop("No probes were found.", call. = FALSE)
        }
        oligos <- rbind(primers, probes)
    } else {
        oligos <- primers
    }
    oligos <- .beautifyOligos(oligos)
    RprimerOligo(oligos)
}

# Helpers ======================================================================

#' Arrange a vector into a matrix of "n-mers"
#'
#' \code{.nmers()} divides a vector into a matrix with \code{n} columns.
#' It is used to generate all possible oligos of a specific length
#' from a DNA sequence.
#'
#' Helper function to \code{.generateOligos()}.
#'
#' @param x A vector.
#'
#' @param n Length of each "mer".
#'
#' @return A matrix with \code{n} columns. The first row will contain
#' \code{x[1:n]}, the second \code{x[2:(n + 1)]}, and so on.
#' In this way, the matrix will contain all possible consecutive
#' sequences of \code{x} of length \code{n}.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' .nmers(c("A", "G", "T", "T", "C", "G"), n = 4)
.nmers <- function(x, n) {
    start <- seq_len(length(x) - n + 1)
    end <- start + n - 1
    nmers <- lapply(start, function(i) x[start[[i]]:end[[i]]])
    do.call("rbind", nmers)
}

#' Count the degeneracy of a DNA sequence
#'
#' \code{.countDegeneracy()} finds the number of unique variants of
#' a DNA sequence with (or without) degenerate bases.
#'
#' Helper function to \code{.generateOligos()}.
#'
#' @param x A DNA sequence (a character vector).
#'
#' @return The number of sequence variants of x.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' .countDegeneracy(c("A", "R", "T", "T", "N", "G"))
.countDegeneracy <- function(x) {
    nNucleotides <- lookup$degeneracy[x]
    prod(nNucleotides)
}

#' Split two matrices, and paste them together
#'
#' Helper function to generate primers according to the "mixed"
#' strategy.
#'
#' @param first
#' The matrix that should appear first.
#'
#' @param last
#' The matrix that should appear last.
#'
#' @param rev
#' If reverse primers should be produced. \code{TRUE} or {FALSE}.
#'
#' @param combine
#' If the two matrices should be pasted together in one. Defaults to \code
#' {TRUE}. Otherwise a list with two matrices will return.
#'
#' @return A matrix.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' .splitAndPaste(t(matrix(rep(1, 10))), t(matrix(rep(2, 10))))
#' .splitAndPaste(t(matrix(rep(1, 10))), t(matrix(rep(2, 10))), combine = FALSE)
#' .splitAndPaste(t(matrix(rep(1, 10))), t(matrix(rep(2, 10))), rev = TRUE)
.splitAndPaste <- function(first, second, rev = FALSE, combine = TRUE) {
    n <- ncol(first)
    small <- seq_len(as.integer(n / 3))
    large <- seq(small[length(small)] + 1, n)
    if (rev) {
        first <- first[, small, drop = FALSE]
        second <- second[, large, drop = FALSE]
    } else {
        first <- first[, seq_along(large), drop = FALSE]
        second <- second[, small + n - length(small), drop = FALSE]
    }
    if (combine) {
        cbind(first, second)
    } else {
        list(first, second)
    }
}

#' Generate oligos of a specific length
#'
#' \code{.generateOligos()} finds all possible oligos of a specific
#' length from an
#' \code{RprimerProfile} object, and returns a list with oligo sequences
#' and additional information.
#'
#' @param x An \code{RprimerProfile} object.
#'
#' @param lengthOligo
#'
#' @return A list.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' .generateOligos(exampleRprimerProfile, lengthOligo = 18)
.generateOligos <- function(x, lengthOligo = 20) {
    oligos <- list()
    oligos$majoritySequence <- .nmers(x$majority, lengthOligo)
    oligos$iupacSequence <- .nmers(x$iupac, lengthOligo)
    oligos$start <- seq_len(nrow(oligos$iupacSequence)) + min(x$position) - 1
    oligos$end <- seq_len(
        nrow(oligos$iupacSequence)
    ) + lengthOligo - 1 + min(x$position) - 1
    oligos$length <- rep(lengthOligo, nrow(oligos$iupacSequence))
    oligos$degeneracy <- apply(oligos$iupacSequence, 1, .countDegeneracy)
    oligos$gapFrequency <- apply(.nmers(x$gaps, lengthOligo), 1, max)
    oligos$coverage <- .nmers(x$coverage, lengthOligo)
    oligos$identity <- .nmers(x$identity, lengthOligo)
    oligos$identityCoverage <- .splitAndPaste(
        oligos$identity, oligos$coverage,
        combine = FALSE
    )
    oligos$coverageIdentity <- .splitAndPaste(
        oligos$coverage, oligos$identity,
        rev = TRUE, combine = FALSE
    )
    oligos$endCoverageFwd <- apply(
        oligos$coverage[
            , (ncol(oligos$coverage) - 5):ncol(oligos$coverage)
        ],
        1, min
    )
    oligos$endCoverageRev <- apply(oligos$coverage[, seq_len(5)], 1, min)
    oligos$identity <- rowMeans(oligos$identity)
    oligos$coverage <- rowMeans(oligos$coverage)
    oligos$method <- rep("ambiguous", nrow(oligos$iupacSequence))
    oligos$roiStart <- rep(
        min(x$position, na.rm = TRUE), nrow(oligos$iupacSequence)
    )
    oligos$roiEnd <- rep(
        max(x$position, na.rm = TRUE), nrow(oligos$iupacSequence)
    )
    oligos
}

#' Mix oligos
#'
#' Helper function to \code{.generateMixedOligos()}
#'
#' @param x An output from \code{.generateOligos()}.
#'
#' @param rev If reverse oligos should be generated.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .generateOligos(exampleRprimerProfile)
#' .mixOligos(x)
.mixOligos <- function(x, rev = FALSE) {
    oligos <- list()
    if (rev) {
        oligos$iupacSequence <- .splitAndPaste(
            x$iupacSequence, x$majoritySequence,
            rev = TRUE
        )
        oligos$coverage <- rowMeans(x$coverageIdentity[[1]])
        oligos$identity <- rowMeans(x$coverageIdentity[[2]])
        oligos$method <- rep("mixedRev", nrow(oligos$iupacSequence))
    } else {
        oligos$iupacSequence <- .splitAndPaste(
            x$majoritySequence, x$iupacSequence
        )
        oligos$coverage <- rowMeans(x$identityCoverage[[2]])
        oligos$identity <- rowMeans(x$identityCoverage[[1]])
        oligos$method <- rep("mixedFwd", nrow(oligos$iupacSequence))
    }
    oligos$degeneracy <- apply(oligos$iupacSequence, 1, .countDegeneracy)
    oligos
}

#' Generate mixed oligos
#'
#' Helper function to \code{.designMixedOligos()}
#'
#' @param x An output from \code{.generateOligos()}.
#'
#' @param rev If reverse oligos should be generated.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .generateOligos(exampleRprimerProfile)
#' .generateMixedOligos(x)
.generateMixedOligos <- function(x, rev = FALSE) {
    oligos <- list()
    oligos <- .mixOligos(x, rev)
    oligos$start <- x$start
    oligos$end <- x$end
    oligos$length <- x$length
    oligos$gapFrequency <- x$gapFrequency
    oligos$endCoverageFwd <- x$endCoverageFwd
    oligos$endCoverageRev <- x$endCoverageRev
    oligos$roiStart <- x$roiStart
    oligos$roiEnd <- x$roiEnd
    order <- c(
        "iupacSequence", "start", "end", "length", "degeneracy", "gapFrequency",
        "identity", "coverage",
        "endCoverageFwd", "endCoverageRev", "method", "roiStart",
        "roiEnd"
    )
    oligos[order]
}

#' Drop items from a list
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .generateOligos(exampleRprimerProfile)
#' .dropItems(x)
.dropItems <- function(x) {
    within(x, rm(
        "majoritySequence", "identityCoverage", "coverageIdentity"
    ))
}

#' Merge two lists with the same structure
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .generateOligos(exampleRprimerProfile)
#' mixedFwd <- .generateMixedOligos(x)
#' mixedRev <- .generateMixedOligos(x, rev = TRUE)
#' .mergeLists(mixedFwd, mixedRev)
.mergeLists <- function(first, second) {
    x <- lapply(names(first), function(y) {
        if (is.matrix(first[[y]])) {
            rbind(first[[y]], second[[y]])
        } else {
            c(first[[y]], second[[y]])
        }
    })
    names(x) <- names(first)
    x
}

#' Design mixed oligos
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .generateOligos(exampleRprimerProfile)
#' .designMixedOligos(x)
.designMixedOligos <- function(x, probe = TRUE) {
    mixedFwd <- .generateMixedOligos(x)
    mixedRev <- .generateMixedOligos(x, rev = TRUE)
    x <- .dropItems(x)
    out <- .mergeLists(mixedFwd, mixedRev)
    if (probe) {
        out <- .mergeLists(out, x)
    }
    out
}

#' Filter oligos upon maximum gap frequency and degeneracy
#'
#' \code{.filterOligos()}
#' removes oligos with gap frequency and degeneracy above the specified
#' thresholds.
#'
#' @param x An output from \code{.generateOligos()}.
#'
#' @inheritParams oligos
#'
#' @return
#' A list with the same structure as \code{x}, but where invalid oligos have
#' been removed.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .generateOligos(exampleRprimerProfile)
#' .filterOligos(x)
.filterOligos <- function(x, maxGapFrequency = 0.1, maxDegeneracy = 4) {
    invalidCharacters <- apply(x$iupacSequence, 1, function(x) {
        any(x == "-") | any(is.na(x))
    })
    invalid <- unique(c(
        which(x$degeneracy > maxDegeneracy),
        which(x$gapFrequency > maxGapFrequency),
        which(invalidCharacters)
    ))
    if (length(invalid > 0)) {
        lapply(x, function(x) {
            if (is.matrix(x)) x[-invalid, , drop = FALSE] else x[-invalid]
        })
    } else {
        x
    }
}

#' Get all variants of a DNA sequence with ambiguous bases
#'
#' Helper function to \code{.getAllVariants()}.
#'
#' @param x A DNA sequence (a character vector).
#'
#' @return All sequence variants of \code{x}. A matrix.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' .expandDegenerates(c("A", "R", "T", "T", "N", "G"))
.expandDegenerates <- function(x) {
    bases <- lapply(x, function(i) {
        degen <- unname(lookup$degenerates[[i]])
        unlist(strsplit(degen, split = ","))
    })
    all <- expand.grid(bases[seq_along(bases)], stringsAsFactors = FALSE)
    all <- as.matrix(all)
    colnames(all) <- NULL
    all
}

#' Turn a list of oligos into a matrix with oligos
#'
#' \code{.makeOligoMatrix()} takes a list of oligos as input,
#' where each element contains a
#' character matrix with all sequence variants of an oligo,
#' and returns a single matrix, where the "oligo ID" of each sequence
#' variant is identified by its rowname.
#'
#' Helper function to \code{.getAllVariants()}.
#'
#' @param x A list with oligos.
#'
#' @return A matrix with all sequence variants, with oligo IDs as rownames.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .filterOligos(.generateOligos(exampleRprimerProfile))
#' x$sequence <- apply(x$iupacSequence, 1, .expandDegenerates)
#' .makeOligoMatrix(x$sequence)
.makeOligoMatrix <- function(x) {
    degeneracy <- vapply(x, nrow, integer(1))
    id <- lapply(seq_along(degeneracy), function(x) rep(x, degeneracy[[x]]))
    id <- unlist(id)
    x <- do.call("rbind", x)
    rownames(x) <- id
    x
}

#' Find the reverse complement of a DNA sequence
#'
#' Helper function to \code{.getAllVariants()}.
#'
#' @param x A matrix with DNA sequence(s).
#'
#' @return The reverse complement of \code{x}, a matrix with the same
#' dimension as \code{x}, and with the same rownames.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' .reverseComplement(matrix(c("A", "R", "T", "T", "N", "G")))
.reverseComplement <- function(x) {
    rc <- x[, rev(seq_len(ncol(x))), drop = FALSE]
    rc[] <- lookup$complement[rc]
    rc
}

#' Identify oligos with a GC-clamp
#'
#' \code{.detectGcClamp()} detects the presence of a GC-clamp.
#' A GC-clamp is identified as two to three G or C:s at
#' the 3'-end (within the last five bases).
#'
#' Helper function to \code{.getAllVariants()}.
#'
#' @param x
#' A numeric matrix, where each row corresponds to a specific oligo.
#' In this matrix, 1 corresponds to G or C, and 0 corresponds to A or T.
#'
#' @param rev
#' If the check should be done in reverse direction.
#'
#' @return
#' A logical vector, where \code{TRUE} indicates the
#' presence of a GC-clamp.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' seq <- matrix(c(1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0))
#' .detectGcClamp(gc, rev = FALSE)
#' .detectGcClamp(gc, rev = TRUE)
.detectGcClamp <- function(x, rev = FALSE) {
    if (rev) {
        end <- x[, seq_len(5), drop = FALSE]
        5 - rowSums(end) >= 2 & 5 - rowSums(end) <= 3 ## Because of complement..
    } else {
        end <- x[, seq((ncol(x) - 4), ncol(x)), drop = FALSE]
        rowSums(end) >= 2 & rowSums(end) <= 3
    }
}

#' Identify oligos with runs of the same nucleotide at the 3' end
#'
#' \code{.detectThreeEndRuns()} checks if the same nucleotide is repeated at
#' at least 3 times at the terminal 3'-end of an oligo (e.g. "AAA")
#'
#' Helper function to \code{.getAllVariants()}.
#'
#' @param x
#' A matrix with DNA sequences.
#'
#' @param rev
#' If the check should be done in reverse direction.
#'
#' @return A logical vector.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' seq <- matrix(c("A", "C", "G", "G", "T", "T", "A", "A"))
#' .detectThreeEndRuns(seq, rev = FALSE)
.detectThreeEndRuns <- function(x, rev = FALSE) {
    if (rev) {
        end <- x[, seq_len(3), drop = FALSE]
    } else {
        end <- x[, seq((ncol(x) - 2), ncol(x)), drop = FALSE]
    }
    apply(end, 1, function(x) {
        all(x == "A") | all(x == "C") | all(x == "T") | all(x == "G")
    })
}

#' Detect mono- and di-nucleotide repeats
#'
#' \code{.detectRepeats()} checks whether an oligo has more than four
#' consecutive "runs" of the same base (e.g. "AAAAA"),
#' or more than three consecutive
#' "runs"  of the same di-nucleotide (e.g. "ATATATAT").
#'
#' Helper function to \code{.getAllVariants()}.
#'
#' @param x A vector with DNA sequences.
#'
#' @return A logical vector, where \code{TRUE} indicates the presence of
#' invalid di- or mononucleotide repeats.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' .detectRepeats(c("ACTTTTCT", "ACTTTTTCT", "ATCTCTCTCA"))
.detectRepeats <- function(x) {
    di <- "(AT){4,}|(TA){4,}|(AC){4,}|(CA){4,}|(AG){4,}|(GA){4,}|(GT){4,}|(TG){4,}|(CG){4,}|(GC){4,}|(CT){4,}|(TC){4,}|)"
    mono <- "([A-Z])\\1\\1\\1\\1"
    vapply(x, function(y) {
        grepl(di, y) | grepl(mono, y)
    }, logical(1))
}

#' Get all variants of oligos with degenerate bases
#'
#' \code{.getAllVariants()} returns all sequence variants of each oligo,
#' both in sense and anti-sense
#' (reverse complement) direction. It calculates GC-content and
#' melting temperature, and provides information on the presence of a GC-clamp,
#' terminal 3'-end runs, mono- and di-nucleotide repeats
#' and terminal five end G:s.
#'
#' Helper function to \code{.designOligos()},
#'
#' @param x An output from \code{.filterOligos()}.
#'
#' @param concPrimer Primer concentration in nM (for tm calculation).
#'
#' @param concProbe Probe concentration in nM (for tm calculation).
#'
#' @param concNa Sodium ion concentration in M (for tm calculation).
#'
#' @inheritParams oligos
#'
#' @return A list.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .filterOligos(.generateOligos(exampleRprimerProfile))
#' .getAllVariants(x)
.getAllVariants <- function(x,
                            concPrimer = 500,
                            concProbe = 250,
                            concNa = 0.05) {
    all <- list()
    all$sequence <- apply(x$iupacSequence, 1, .expandDegenerates)
    ## If there is only one variant of each oligo,
    ## (and apply returns a matrix instead of a list):
    if (is.matrix(all$sequence)) {
        all$sequence <- t(all$sequence)
        all$sequence <- lapply(seq_len(nrow(all$sequence)), function(i) {
            all$sequence[i, , drop = FALSE]
        })
    }
    all$sequence <- .makeOligoMatrix(all$sequence)
    all$sequenceRc <- .reverseComplement(all$sequence)
    gc <- all$sequence == "C" | all$sequence == "G"
    n <- rowSums(
        all$sequence == "A" | all$sequence == "C" |
            all$sequence == "G" | all$sequence == "T"
    )
    all$gcContent <- rowSums(gc) / n
    all$gcClampFwd <- .detectGcClamp(gc)
    all$gcClampRev <- .detectGcClamp(gc, rev = TRUE)
    all$threeEndRunsFwd <- .detectThreeEndRuns(all$sequence)
    all$threeEndRunsRev <- .detectThreeEndRuns(all$sequence, rev = TRUE)
    all$fiveEndGPlus <- all$sequence[, 1] == "G"
    all$fiveEndGMinus <- all$sequence[, ncol(all$sequence)] == "C"
    tmParam <- .tmParameters(all$sequence, concNa)
    all$tmPrimer <- apply(tmParam, 1, function(x) .tm(x, concPrimer))
    all$tmProbe <- apply(tmParam, 1, function(x) .tm(x, concProbe))
    all$sequence <- apply(all$sequence, 1, paste, collapse = "")
    all$repeats <- .detectRepeats(all$sequence)
    all$sequenceRc <- apply(all$sequenceRc, 1, paste, collapse = "")
    lapply(all, function(x) unname(split(unname(x), f = as.integer(names(x)))))
}

#' Calculate mean values and ranges for GC-content and Tm
#'
#' When all sequence variants of each oligo are generated with
#' \code{.getAllVariants()}, the next step is to compute mean values and ranges
#' for GC-content and melting temperature.
#'
#' Helper function to \code{.designOligos()}.
#'
#' @param x An output from \code{.getAllVariants()}.
#'
#' @return A data frame with mean values and ranges.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .getAllVariants(.filterOligos(.generateOligos(exampleRprimerProfile)))
#' .getMeanAndRange(x)
.getMeanAndRange <- function(x) {
    x <- x[c("gcContent", "tmPrimer", "tmProbe")]
    means <- lapply(x, function(y) {
        vapply(y, function(z) {
            sum(z) / length(z)
        }, double(1))
    })
    means <- do.call("cbind.data.frame", means)
    names(means) <- paste0(names(means), "Mean")
    ranges <- lapply(x, function(y) {
        vapply(y, function(z) {
            max(z) - min(z)
        }, double(1))
    })
    ranges <- do.call("cbind.data.frame", ranges)
    names(ranges) <- paste0(names(ranges), "Range")
    cbind(means, ranges)
}

#' Turn a list with oligos into a data frame
#'
#' Helper function to \code{.designOligos()}.
#'
#' @param x A list with oligos.
#'
#' @return A data frame with oligos.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .dropItems(.filterOligos(.generateOligos(exampleRprimerProfile)))
#' .makeOligoDf(x)
.makeOligoDf <- function(x) {
    x <- within(x, rm("gapFrequency"))
    x$iupacSequenceRc <- .reverseComplement(x$iupacSequence)
    x$iupacSequence <- apply(x$iupacSequence, 1, paste, collapse = "")
    x$iupacSequenceRc <- apply(x$iupacSequenceRc, 1, paste, collapse = "")
    do.call("cbind.data.frame", x)
}

#' Design oligos
#'
#' Helper function to \code{oligos()}.
#'
#' @inheritParams oligos
#'
#' @return A data frame with oligos.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' .designOligos(exampleRprimerProfile)
.designOligos <- function(x,
                          lengthOligo = 18:22,
                          maxGapFrequency = 0.1,
                          maxDegeneracy = 4,
                          concPrimer = 500,
                          designStrategyPrimer = "ambiguous",
                          probe = TRUE,
                          concProbe = 250,
                          concNa = 0.05) {
    allOligos <- lapply(lengthOligo, function(i) {
        iupacOligos <- .generateOligos(x, lengthOligo = i)
        if (designStrategyPrimer == "mixed") {
            iupacOligos <- .designMixedOligos(iupacOligos, probe)
        } else {
            iupacOligos <- .dropItems(iupacOligos)
        }
        iupacOligos <- .filterOligos(
            iupacOligos, maxGapFrequency, maxDegeneracy
        )
        nOligos <- vapply(iupacOligos, length, integer(1))
        if (all(nOligos == 0L)) {
            stop("No primers were found.", call. = FALSE)
        }
        allVariants <- .getAllVariants(
            iupacOligos, concPrimer, concProbe, concNa
        )
        meansAndRanges <- .getMeanAndRange(allVariants)
        allVariants <- data.frame(do.call("cbind", allVariants))
        iupacOligos <- .makeOligoDf(iupacOligos)
        cbind(iupacOligos, meansAndRanges, allVariants)
    })
    do.call("rbind", allOligos)
}

#' Check if vectors in a list are within a specified range
#'
#' Helper function to \code{.filterPrimers()} and \code{.filterProbes()}.
#'
#' @param x A list.
#'
#' @param range The specified range. A numeric vector of length two.
#'
#' @return A list with logical vectors.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .getAllVariants(.filterOligos(.generateOligos(exampleRprimerProfile)))
#' .isWithinRange(x$gcContent, c(0.4, 0.6))
.isWithinRange <- function(x, range) {
    lapply(x, function(y) y >= min(range) & y <= max(range))
}

#' Convert a data frame with lists to a list of matrices
#'
#' Helper function to \code{.checkAllPrimerVariants()} and
#' \code{.checkAllProbeVariants()}
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .designOligos(exampleRprimerProfile)
#' gcInRange <- .isWithinRange(x$gcContent, c(0.4, 0.6))
#' x <- cbind(x, data.frame(cbind(gcInRange)))
#' .convertToMatrices(x["gcInRange"])
.convertToMatrices <- function(x) {
    lapply(seq_len(nrow(x)), function(i) {
        y <- lapply(x[i, , drop = FALSE], unlist)
        do.call("cbind", y)
    })
}

#' Check if design constraints are fulfilled
#'
#' \code{.checkValidity()} takes a list with logical matrices as input.
#' Within these matrices, each row represents a unique sequence variant
#' and each column
#' represents a design constraint (e.g. the presence of GC-clamp, or
#' dinucleotide repeats, etc.). It inverts columns with "undesired" criteria
#' (e.g. the presence of dinucleotide repeats), computes row and column means,
#' and checks
#' whether the row and column means are equal to or higher than
#' a specified acceptance threshold.
#'
#' Column means represent the proportion
#' of sequence variants that fulfill a specific criteria (e.g. GC-clamp),
#' and row means represent the proportion of the desired design criteria that
#' are fulfilled by specific sequence variants.
#'
#' Helper function to \code{.filterPrimers()} and \code{filterProbes()}.
#'
#' @param x A list with logical matrices.
#'
#' @param rowThreshold
#' Minimum accepted value for row means.
#'
#' @param colThreshold Minimum accepted value for column means.
#'
#' @return A logical vector.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .designOligos(exampleRprimerProfile)
#' gcInRange <- .isWithinRange(x$gcContent, c(0.4, 0.6))
#' x <- cbind(x, data.frame(cbind(gcInRange)))
#' check <- .convertToMatrices(x[c("gcInRange", "repeats", "threeEndRunsFwd")])
#' .isValid(check, rowThreshold = 0.5, colThreshold = 0.5)
.isValid <- function(x, rowThreshold, colThreshold) {
    valid <- vapply(x, function(y) {
        toInvert <- c(
            "repeats", "threeEndRunsFwd", "threeEndRunsRev",
            "fiveEndGPlus", "fiveEndGMinus"
        )
        select <- colnames(y) %in% toInvert
        y[, select] <- !y[, select]
        col <- colMeans(y)
        row <- rowMeans(y)
        all(col >= colThreshold) & all(row >= rowThreshold)
    }, logical(1))
    valid
}

#' Check all primer variants
#'
#' Helper function to \code{.filterPrimers()}
#'
#' @param x A data frame (see examples)
#'
#' @return A data frame
#'
#' @inheritParams oligos
#'
#' @param rowThreshold
#' Minimum proportion of the specified design constraints that must be
#' met by each sequence variant.
#'
#' @param colThreshold
#' Minimum proportion of the sequence variants that must meet each
#' design constraint.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .designOligos(exampleRprimerProfile)
#' gcInRange <- .isWithinRange(x$gcContent, c(0.4, 0.6))
#' tmInRange <- .isWithinRange(x$tmPrimer, c(55, 65))
#' x <- cbind(x, data.frame(cbind(tmInRange, gcInRange)))
#' .checkAllPrimerVariants(x)
.checkAllPrimerVariants <- function(x,
                                    gcClampPrimer = TRUE,
                                    avoidThreeEndRunsPrimer = TRUE,
                                    rowThreshold = 1,
                                    colThreshold = 1) {
    selectFwd <- c("repeats", "tmInRange", "gcInRange")
    selectRev <- c("repeats", "tmInRange", "gcInRange")
    if (gcClampPrimer) {
        selectFwd <- c(selectFwd, "gcClampFwd")
        selectRev <- c(selectRev, "gcClampRev")
    }
    if (avoidThreeEndRunsPrimer) {
        selectFwd <- c(selectFwd, "threeEndRunsFwd")
        selectRev <- c(selectRev, "threeEndRunsRev")
    }
    xFwd <- x[selectFwd]
    xFwd <- .convertToMatrices(xFwd)
    okFwd <- .isValid(xFwd, rowThreshold, colThreshold)
    xRev <- x[selectRev]
    xRev <- .convertToMatrices(xRev)
    okRev <- .isValid(xRev, rowThreshold, colThreshold)
    x <- cbind(x, okFwd, okRev)
    x[x$okFwd | x$okRev, , drop = FALSE]
}

#' Find oligos that pass the criteria for being a primer
#'
#' Helper function to \code{oligos()}
#'
#' @param x An output from \code{.designOligos()}
#'
#' @inheritParams oligos
#'
#' @param rowThreshold
#' Minimum proportion of the specified design constraints that must be
#' met by each sequence variant.
#'
#' @param colThreshold
#' Minimum proportion of the sequence variants that must meet each
#' design constraint.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .designOligos(exampleRprimerProfile)
#' .filterPrimers(x)
.filterPrimers <- function(x,
                           lengthPrimer = 18:22,
                           maxDegeneracyPrimer = 4,
                           gcClampPrimer = TRUE,
                           avoidThreeEndRunsPrimer = TRUE,
                           minThreeEndCoveragePrimer = 0.98,
                           gcRangePrimer = c(0.45, 0.55),
                           tmRangePrimer = c(55, 65),
                           designStrategyPrimer = "ambiguous",
                           colThreshold = 0.75,
                           rowThreshold = 0.75) {
    x <- x[x$length %in% lengthPrimer, , drop = FALSE]
    x <- x[x$degeneracy <= maxDegeneracyPrimer, , drop = FALSE]
    x <- x[
        x$endCoverageFwd >= minThreeEndCoveragePrimer |
            x$endCoverageFwd >= minThreeEndCoveragePrimer, ,
        drop = FALSE
    ]
    gcInRange <- .isWithinRange(x$gcContent, gcRangePrimer)
    tmInRange <- .isWithinRange(x$tmPrimer, tmRangePrimer)
    x <- cbind(x, data.frame(cbind(tmInRange, gcInRange)))
    x <- .checkAllPrimerVariants(
        x,
        gcClampPrimer,
        avoidThreeEndRunsPrimer,
        rowThreshold,
        colThreshold
    )
    fwd <- x$endCoverageFwd >= minThreeEndCoveragePrimer & x$okFwd
    rev <- x$endCoverageRev >= minThreeEndCoveragePrimer & x$okRev
    x <- cbind(x, fwd, rev)
    x$rev[x$method == "mixedFwd" & x$rev] <- FALSE
    x$fwd[x$method == "mixedRev" & x$fwd] <- FALSE
    x <- x[x$fwd | x$rev, , drop = FALSE]
    remove <- c(
        "gcInRange", "tmInRange", "endCoverageFwd", "endCoverageRev",
        "okFwd", "okRev", "tmProbeMean", "tmProbeRange", "tmProbe"
    )
    x <- x[!names(x) %in% remove]
    oldnames <- c("tmPrimerMean", "tmPrimerRange", "tmPrimer")
    newnames <- c("tmMean", "tmRange", "tm")
    names(x)[names(x) %in% oldnames] <- newnames
    type <- rep("primer", nrow(x))
    cbind(type, x)
}

#' Check all probe variants
#'
#' Helper function to \code{.filterPrimers()}
#'
#' @param x
#'
#' @inheritParams oligos
#'
#' @param rowThreshold
#' Minimum proportion of the specified design constraints that must be
#' met by each sequence variant.
#'
#' @param colThreshold
#' Minimum proportion of the sequence variants that must meet each
#' design constraint.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .designOligos(exampleRprimerProfile)
#' gcInRange <- .isWithinRange(x$gcContent, c(0.4, 0.6))
#' tmInRange <- .isWithinRange(x$tmPrimer, c(55, 65))
#' x <- cbind(x, data.frame(cbind(tmInRange, gcInRange)))
#' .checkAllProbeVariants(x)
.checkAllProbeVariants <- function(x,
                                   avoidFiveEndGProbe,
                                   rowThreshold,
                                   colThreshold) {
    selectFwd <- c("repeats", "tmInRange", "gcInRange")
    selectRev <- c("repeats", "tmInRange", "gcInRange")
    if (avoidFiveEndGProbe) {
        selectFwd <- c(selectFwd, "fiveEndGPlus")
        selectRev <- c(selectRev, "fiveEndGMinus")
    }
    xFwd <- x[selectFwd]
    xFwd <- .convertToMatrices(xFwd)
    fwd <- .isValid(xFwd, rowThreshold, colThreshold)
    xRev <- x[selectRev]
    xRev <- .convertToMatrices(xRev)
    rev <- .isValid(xRev, rowThreshold, colThreshold)
    x <- cbind(x, fwd, rev)
    x[x$fwd | x$rev, , drop = FALSE]
}

#' Find oligos that pass the criteria for being a probe
#'
#' Helper function to \code{oligos()}
#'
#' @param x An output from \code{.generateOligos()}
#'
#' @inheritParams oligos
#'
#' @param rowThreshold
#' Minimum proportion of the specified design constraints that must be
#' met by each sequence variant.
#'
#' @param colThreshold
#' Minimum proportion of the sequence variants that must meet each
#' design constraint.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .designOligos(exampleRprimerProfile)
#' .filterProbes(x)
.filterProbes <- function(x,
                          lengthProbe = 18:22,
                          maxDegeneracyProbe = 4,
                          avoidFiveEndGProbe = TRUE,
                          gcRangeProbe = c(0.45, 0.55),
                          tmRangeProbe = c(55, 65),
                          rowThreshold = 0.75,
                          colThreshold = 0.75) {
    x <- x[x$method == "ambiguous", , drop = FALSE]
    x <- x[x$length %in% lengthProbe, , drop = FALSE]
    x <- x[x$degeneracy <= maxDegeneracyProbe, , drop = FALSE]
    gcInRange <- .isWithinRange(x$gcContent, gcRangeProbe)
    tmInRange <- .isWithinRange(x$tmProbe, tmRangeProbe)
    x <- cbind(x, data.frame(cbind(tmInRange, gcInRange)))
    x <- .checkAllProbeVariants(
        x,
        avoidFiveEndGProbe,
        rowThreshold,
        colThreshold
    )
    remove <- c(
        "gcInRange", "tmInRange", "tmPrimerMean",
        "endCoverageFwd", "endCoverageRev", "tmPrimerRange", "tmPrimer"
    )
    x <- x[!names(x) %in% remove]
    oldnames <- c("tmProbeMean", "tmProbeRange", "tmProbe")
    newnames <- c("tmMean", "tmRange", "tm")
    names(x)[names(x) %in% oldnames] <- newnames
    type <- rep("probe", nrow(x))
    cbind(type, x)
}

#' Beautify oligo data
#'
#' \code{.beautifyOligos()} drops unnecessary columns and sorts oligos based
#' on their start position.
#'
#' Helper function to \code{oligos()}.
#'
#' @param x A data frame with oligos.
#'
#' @return A data frame with oligos.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- .filterPrimers(.designOligos(exampleRprimerProfile))
#' .beautifyOligos(x)
.beautifyOligos <- function(x) {
    keep <- c(
        "type", "fwd", "rev", "start", "end", "length",
        "iupacSequence", "iupacSequenceRc", "identity",
        "coverage", "degeneracy", "gcContentMean", "gcContentRange",
        "tmMean", "tmRange", "sequence",
        "sequenceRc", "gcContent", "tm", "method", "roiStart",
        "roiEnd"
    )
    x <- x[keep]
    x <- x[order(x$start), ]
    rownames(x) <- NULL
    x
}
