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
#' Primer length. A numeric vector [14, 30],
#' defaults to \code{18:22}.
#'
#' @param maxDegeneracyPrimer
#' Maximum number of variants of each primer. A number [1, 32], defaults to 4.
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
#' Seven bases ######################################################################
#'
#' @param gcRangePrimer
#' GC-content range for primers (proportion, not percent).
#' A numeric vector [0, 1], defaults to \code{c(0.40, 0.60)}.
#'
#' @param tmRangePrimer
#' Tm range for primers.
#' A numeric vector [30, 90], defaults to \code{c(55, 65)}.
#'
#' @param concPrimer
#' Primer concentration in nM, for Tm calculation. A number
#' [20, 2000], defaults to 500.
#'
#' @param probe
#' If probes should be designed. \code{TRUE} or \code{FALSE},
#' defaults to \code{TRUE}.
#'
#' @param lengthProbe
#' Probe length. A numeric vector [14, 30],
#' defaults to \code{18:22}.
#'
#' @param maxDegeneracyProbe
#' Maximum number of variants of each probe. A number [1, 32], defaults to 4.
#'
#' @param avoidFiveEndGProbe
#' If probes with G
#' at the 5'-end should be avoided. \code{TRUE} or \code{FALSE},
#' defaults to \code{TRUE}.
#'
#' @param gcRangeProbe
#' GC-content range for probes (proportion, not %). A numeric vector [0, 1],
#' defaults to \code{c(0.40, 0.60)}.
#'
#' @param tmRangeProbe
#' Tm range for probes.
#' A numeric vector [30, 90], defaults to \code{c(55, 65)}.
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
#'   \item{iupacSequence}{Oligo sequence, with wobble bases (if any).}
#'   \item{iupaSequenceRc}{The reverse complement of the iupacSequence.}
#'   \item{identity}{Average identity score of the oligo, can range from 0 to 1.
#'     The identity is the proportion of the most common base at each position
#'     in the input alignment.}
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
#'   \item{tm}{tm of all sequence variants.}
#'   \item{roiStart}{First position of the input \code{RprimerProfile} object
#'     (roi = region of interest).}
#'   \item{roiEnd}{Last position of the input \code{RprimerProfile} object.}
#' }
#'
#' @section Oligos with low sequence complexity:
#'
#' Oligos with more than four consecutive runs of the same
#' nucleotide (e.g. "AAAAA") and/or more than three consecutive runs
#' of the same di-nucleotide (e.g. "TATATATA") are considered invalid.
#' This check is done on all sequence variants of each oligo. ################################
#'
#' @section Oligos with high degeneracy:
#'
#' Each sequence variant of a degenerate oligo must fulfill at least
#' 83 % of the specified design criteria,
#' and each design criteria must be fulfilled by at least 5/6 of
#' the sequence variants. Otherwise the oligo will be considered as invalid. ###############
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
#'        maxDegeneracyPrimer = 32,
#'        probe = FALSE
#' )
oligos <- function(x,
                   maxGapFrequency = 0.05,
                   lengthPrimer = 18:22,
                   maxDegeneracyPrimer = 4,
                   gcClampPrimer = TRUE,
                   avoidThreeEndRunsPrimer = TRUE,
                   minThreeEndCoveragePrimer = 0.98,
                   gcRangePrimer = c(0.40, 0.60),
                   tmRangePrimer = c(55, 65),
                   concPrimer = 500,
                   probe = TRUE,
                   lengthProbe = 18:22,
                   maxDegeneracyProbe = 4,
                   avoidFiveEndGProbe = TRUE,
                   gcRangeProbe = c(0.40, 0.60),
                   tmRangeProbe = c(55, 70),
                   concProbe = 250,
                   concNa = 0.05) {
    if (!methods::is(x, "RprimerProfile")) {
        stop("'x' must be an RprimerProfile object.", call. = FALSE)
    }
    if (!(maxGapFrequency >= 0 && maxGapFrequency <= 1)) {
        stop("'lengthPrimer' must be from 0 to 1.", call. = FALSE)
    }
    if (!(min(lengthPrimer) >= 14 && max(lengthPrimer) <= 30)) {
        stop("'lengthPrimer' must be from 14 to 30.", call. = FALSE)
    }
    if (!(maxDegeneracyPrimer >= 1 && maxDegeneracyPrimer <= 32)) {
        stop("'maxDegeneracyPrimer' must be from 1 to 32.", call. = FALSE)
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
        stop("'concPrimer' must be from 20 nM to 2000 nM.", call. = FALSE)
    }
    if (!is.logical(probe)) {
        stop("'probe' must be set to TRUE or FALSE", call. = FALSE)
    }
    if (!(min(lengthProbe) >= 14 && max(lengthProbe) <= 30)) {
        stop("'lengthProbe' must be from 14 to 30.", call. = FALSE)
    }
    if (!(maxDegeneracyProbe >= 1 && maxDegeneracyProbe <= 32)) {
        stop("'maxDegeneracyProbe' must be from 1 to 32.", call. = FALSE)
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
        stop("'concProbe' must be from 20 nM to 2000 nM.", call. = FALSE)
    }
    if (!(concNa >= 0.01 && concNa <= 1)) {
        stop("'concNa' must be from 0.01 to 1 M.", call. = FALSE)
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
        rowThreshold = 1,
        colThreshold = 1
    )
    if (nrow(primers) == 0) {
        stop("No primers were found.", call. = FALSE)
    }
    if (probe) {
        probes <- .filterProbes(oligos,
            lengthProbe,
            maxDegeneracyProbe,
            avoidFiveEndGProbe,
            gcRangeProbe,
            tmRangeProbe,
            rowThreshold = 0.75,
            colThreshold = 0.75
        )
        if (nrow(probes) == 0) {
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
#' @param n Length of each "mer" (a positive integer).
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
#' @return The number of sequence variants of x (an integer).
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

#' Generate oligos of a specific length
#'
#' \code{.generateOligos()} is the first step of the oligo-design process.
#' It finds all possible oligos of a specific length from an
#' \code{RprimerProfile} object, and returns a list containing
#' start and end position,
#' length, IUPAC sequence (the oligo DNA sequence with wobble bases), degeneracy
#' (number of variants of each oligo), maximum gap frequency,
#' mean overall identity, and minimum 3'-end identity at
#' both forward and reverse direction.
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
    oligos$iupacSequence <- .nmers(x$iupac, lengthOligo)
    oligos$start <- seq_len(nrow(oligos$iupacSequence)) + min(x$position) - 1
    oligos$end <- seq_len(
        nrow(oligos$iupacSequence)
    ) + lengthOligo - 1 + min(x$position) - 1
    oligos$length <- rep(lengthOligo, nrow(oligos$iupacSequence))
    oligos$degeneracy <- apply(oligos$iupacSequence, 1, .countDegeneracy)
    oligos$gapFrequency <- apply(.nmers(x$gaps, lengthOligo), 1, max)
    oligos$coverage <- .nmers(x$coverage, lengthOligo)
    oligos$endCoverageFwd <- apply(
        oligos$coverage[
            , (ncol(oligos$coverage) - 6):ncol(oligos$coverage)],
        1, min
    )
    oligos$endCoverageRev <- apply(oligos$coverage[, seq_len(6)], 1, min)
    oligos$coverage <- rowMeans(oligos$coverage)
    oligos$roiStart <- rep(
        min(x$position, na.rm = TRUE), nrow(oligos$iupacSequence)
    )
    oligos$roiEnd <- rep(
        max(x$position, na.rm = TRUE), nrow(oligos$iupacSequence)
    )
    oligos
}

#' Remove oligos with too high gap frequency and degeneracy
#'
#' \code{.filterOligos()} is the second step of the oligo-design process. It
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
.filterOligos <- function(x, maxGapFrequency = 0.1, maxDegeneracy = 4) {
    invalid <- unique(c(
        which(x$degeneracy > maxDegeneracy),
        which(x$gapFrequency > maxGapFrequency)
    ))
    if (length(invalid > 0)) {
        lapply(x, function(x) {
            if (is.matrix(x)) x[-invalid, , drop = FALSE] else x[-invalid]
        })
    } else {
        x
    }
}

#' Get all variants of a DNA sequence with wobble bases
#'
#' Helper function to \code{.getAllVariants()}.
#'
#' @param x A DNA sequence, i.e. a character vector.
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
#' \code{.makeOligoMatrix()} is part of a "workaround" to avoid loops when
#' finding e.g. the reverse complement, the presence of GC-clamp etc.
#' of many oligos.
#'
#' It takes a list of DNA sequences as input, where each element contains a
#' character matrix with all sequence variants of a specific degenerate oligo,
#' and returns a single matrix, where the "oligo-belonging" of each sequence
#' is identified by a rowname ID.
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
#' \code{.detectGcClamp()} detects the presence of a GC-clamp on oligos
#' (good to have on primers).
#'
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
#' A logical vector of length \code{nrow(x)}, where \code{TRUE} indicates the
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
        ifelse(5 - rowSums(end) >= 2 & 5 - rowSums(end) <= 3, TRUE, FALSE)
    } else {
        end <- x[, seq((ncol(x) - 4), ncol(x)), drop = FALSE]
        ifelse(rowSums(end) >= 2 & rowSums(end) <= 3, TRUE, FALSE)
    }
}

#' Identify oligos with runs of the same nucleotide at the 3' end
#'
#' \code{.detectThreeEndRuns()} detects if the same nucleotide is repeated at
#' at least 3 times at the terminal 3'-end of an oligo (e.g. "AAA")
#' (bad to have on primers).
#'
#' Helper function to \code{.getAllVariants()}.
#'
#' @param x
#' A matrix with DNA sequences.
#'
#' @param rev
#' If the check should be done in reverse direction.
#'
#' @return
#' A logical vector of length \code{nrow(x)}, where \code{TRUE} indicates the
#' presence of a 3'-end run.
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
    vapply(x, function(i) {
        ifelse(grepl(di, i) | grepl(mono, i), TRUE, FALSE)
    }, logical(1))
}

#' Get all variants of oligos with degenerate bases
#'
#' \code{.getAllVariants()} is the third step of the oligo-design process.
#' It returns all sequence variants of each oligo, both in sense and anti-sense
#' (reverse complement) direction. It also calculates GC-content and
#' melting temperature, and provides information on the presence of a GC-clamp
#' (good for primers),
#' terminal 3'-end runs (bad for primers), mono- and di-nucleotide repeats
#' (bad for primers and probes) and terminal five end G:s (bad for probes).
#'
#' Helper function to \code{.designOligos()},
#'
#' @param x An output from \code{.filterOligos()}.
#'
#' @inheritParams oligos
#'
#' @return
#' A list.
#'
#' @keywords internal
#'
#' @noRd
.getAllVariants <- function(x,
                            concPrimer = 500,
                            concProbe = 250,
                            concNa = 0.05) {
    all <- list()
    all$sequence <- apply(x$iupacSequence, 1, .expandDegenerates)
    ## some kind of workaround if there is only one variant of each oligo,
    ## and apply returns a matrix instead of a list...
    if (!is.list(all$sequence)) {
        all$sequence <- t(all$sequence)
        all$sequence <- lapply(seq_len(nrow(all$sequence)), function(i) {
            all$sequence[i, , drop = FALSE]
        })
    }
    all$sequence <- .makeOligoMatrix(all$sequence)
    all$sequenceRc <- .reverseComplement(all$sequence)
    gc <- ifelse(all$sequence == "C" | all$sequence == "G", 1, 0)
    n <- rowSums(ifelse(
        all$sequence == "A" | all$sequence == "C" |
            all$sequence == "G" | all$sequence == "T", 1, 0
    ))
    all$gcContent <- rowSums(gc) / n
    all$gcClampFwd <- .detectGcClamp(gc)
    all$gcClampRev <- .detectGcClamp(gc, rev = TRUE)
    all$threeEndRunsFwd <- .detectThreeEndRuns(all$sequence)
    all$threeEndRunsRev <- .detectThreeEndRuns(all$sequence, rev = TRUE)
    all$fiveEndGPlus <- ifelse(all$sequence[, 1] == "G", TRUE, FALSE)
    all$fiveEndGMinus <- ifelse(
        all$sequence[, ncol(all$sequence)] == "C", TRUE, FALSE
    )
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
#' @param x An output from \code{.generateOligos()} or \code{.filterOligos()}.
#'
#' @return A data frame with oligos.
#'
#' @keywords internal
#'
#' @noRd
.makeOligoDf <- function(x) {
    gapFrequency <- NULL ## Just to avoid cmd check note
    x <- within(x, rm(gapFrequency))
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
.designOligos <- function(x,
                          lengthOligo = 18:22,
                          maxGapFrequency = 0.1,
                          maxDegeneracy = 4,
                          concPrimer = 500,
                          concProbe = 250,
                          concNa = 0.05) {
    allOligos <- lapply(lengthOligo, function(i) {
        iupacOligos <- .generateOligos(x, lengthOligo = i)
        iupacOligos <- .filterOligos(iupacOligos,
            maxGapFrequency = maxGapFrequency,
            maxDegeneracy = maxDegeneracy
        )
        nOligos <- vapply(iupacOligos, length, integer(1))
        if (all(nOligos == 0)) {
            stop("No primers were found.", call. = FALSE)
        }
        allVariants <- .getAllVariants(
            iupacOligos,
            concPrimer,
            concProbe,
            concNa
        )
        meansAndRanges <- .getMeanAndRange(allVariants)
        allVariants <- data.frame(do.call("cbind", allVariants))
        iupacOligos <- .makeOligoDf(iupacOligos)
        cbind(iupacOligos, meansAndRanges, allVariants)
    })
    do.call("rbind", allOligos)
}

#' Check if vectors in a list are within a specificed range
#'
#' Helper function to \code{.filterPrimers()} and \code{.filterProbes()}.
#'
#' @param x A list.
#'
#' @param range The specified range. A numeric vector of length two.
#'
#' @keywords internal
#'
#' @noRd
.isWithinRange <- function(x, range) {
    lapply(x, function(y) {
        ifelse(y >= min(range) & y <= max(range), TRUE, FALSE)
    })
}

#' Convert a data frame with lists to a list of matrices
#'
#' Helper function to \code{.filterPrimers()} and \code{.filterProbes()}.
#'
#' @keywords internal
#'
#' @noRd
.convertToMatrices <- function(x) {
    lapply(seq_len(nrow(x)), function(i) {
        y <- lapply(x[i, ], unlist)
        y <- do.call("cbind", y)
        if (!is.matrix(y)) t(matrix(y)) else y # if one row, test this  ...
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
.isValid <- function(x, rowThreshold, colThreshold) {
    valid <- vapply(x, function(y) {
        toInvert <- c(
            "repeats", "threeEndRunsFwd", "threeEndRunsRev",
            "fiveEndGPlus", "fiveEndGMinus"
        )
        select <- colnames(y) %in% toInvert
        y[, select] <- as.logical(1 - y[, select])
        col <- colMeans(y)
        row <- rowMeans(y)
        if (all(col >= colThreshold) & all(row >= rowThreshold)) TRUE else FALSE
    }, logical(1))
    valid
}

#' Check all primer variants
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
.checkAllPrimerVariants <- function(x,
                                    gcClampPrimer,
                                    avoidThreeEndRunsPrimer,
                                    rowThreshold,
                                    colThreshold) {
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
.filterPrimers <- function(x,
                           lengthPrimer = 18:22,
                           maxDegeneracyPrimer = 4,
                           gcClampPrimer = TRUE,
                           avoidThreeEndRunsPrimer = TRUE,
                           minThreeEndCoveragePrimer = 0.98,
                           gcRangePrimer = c(0.45, 0.55),
                           tmRangePrimer = c(55, 65),
                           colThreshold = 0.75,
                           rowThreshold = 0.75) {
    x <- x[
        x$length >= min(lengthPrimer) & x$length <= max(lengthPrimer), , ##
        drop = FALSE
    ]
    x <- x[x$degeneracy <= maxDegeneracyPrimer, , drop = FALSE]
    if (minThreeEndCoveragePrimer) {
        x <- x[
            x$endCoverageFwd >= minThreeEndCoveragePrimer |
                x$endCoverageFwd >= minThreeEndCoveragePrimer, , drop = FALSE
        ]
    }
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
    if (minThreeEndCoveragePrimer) {
        fwd <- x$endCoverageFwd >= minThreeEndCoveragePrimer & x$okFwd
        rev <- x$endCoverageRev >= minThreeEndCoveragePrimer & x$okRev
    } else {
        fwd <- x$okFwd
        rev <- x$okRev
    }
    x <- cbind(x, fwd, rev)
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
.filterProbes <- function(x,
                          lengthProbe = 18:22,
                          maxDegeneracyProbe = 4,
                          avoidFiveEndGProbe = TRUE,
                          gcRangeProbe = c(0.45, 0.55),
                          tmRangeProbe = c(55, 65),
                          rowThreshold = 0.75,
                          colThreshold = 0.75) {
    x <- x[
        x$length >= min(lengthProbe) & x$length <= max(lengthProbe), ,
        drop = FALSE
    ]
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
.beautifyOligos <- function(x) {
    keep <- c(
        "type", "fwd", "rev", "start", "end", "length",
        "iupacSequence", "iupacSequenceRc",
        "coverage", "degeneracy", "gcContentMean", "gcContentRange",
        "tmMean", "tmRange", "sequence",
        "sequenceRc", "gcContent", "tm", "roiStart",
        "roiEnd"
    )
    x <- x[keep]
    x <- x[order(x$start), ]
    rownames(x) <- NULL
    x
}
