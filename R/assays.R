#' Design (RT)-PCR assays
#'
#' \code{assays()} combines forward, reverse primers
#' and probes to (RT)-PCR assays from an \code{RprimerOligo} object.
#'
#' @param x An \code{RprimerOligo} object, with or without probes.
#'
#' @param lengthRange
#' Amplicon length range, a numeric vector [40, 5000], defaults to
#' \code{c(65, 120)}.
#'
#' @param tmDiffPrimers
#' Maximum allowed difference between the mean tm of the forward primer
#' and the mean tm of the reverse primer
#' (absolute value). A number [0, 20], defaults to 5.
#'
#' @param tmDiffPrimersProbe
#' Acceptable difference between the mean tm of the probe and the mean tm
#' of the primers. A numeric vector [-20, 20],
#' defaults to \code{c(-10, 10)}.
#' A negative tm difference
#' means that the mean tm of the probe is lower than the mean tm of the
#' primers.
#'
#' @return
#' An \code{RprimerAssay} object.
#' An error message will return if no assays are found.
#'
#' The object contains the following information:
#'
#' \describe{
#'   \item{start}{Position where the assay starts.}
#'   \item{end}{Position where the assay ends.}
#'   \item{ampliconLength}{Length of the amplicon.}
#'   \item{tmDifferencePrimer}{Difference between
#'   the mean tm of the forward primer and the mean tm of the reverse primer,
#'   absolute value.}
#'   \item{deltaGDifferencePrimer}{Difference between
#'   the mean delta G of the forward primer and the mean delta G of
#'   the reverse primer, absolute value.}
#'   \item{totalDegeneracy}{Total number of oligos in the assay.}
#'   \item{score}{Summarized oligo score. The lowest, and best,
#'   possible score is 0. The highest possible score is 24 for assays
#'   with only primers,
#'   and 36 for assays with probes.
#'   See \code{?oligos} for more information about the scoring system.}
#'   \item{startFwd}{Position where the forward primer starts.}
#'   \item{endFwd}{Position where the forward primer ends.}
#'   \item{lengthFwd}{Length of the forward primer.}
#'   \item{iupacSequenceFwd}{IUPAC sequence of the forward primer.}
#'   \item{coverageFwd}{Average coverage of the forward primer.}
#'   \item{degeneracyFwd}{Number of variants of the forward primer.}
#'   \item{gcContentMeanFwd}{Mean GC-content of the forward primer.}
#'   \item{gcContentRangeFwd}{Range in GC-content of the forward primer.}
#'   \item{tmMeanFwd}{Mean tm of the forward primer.}
#'   \item{tmRangeFwd}{Range in tm of the forward primer.}
#'   \item{deltaGMeanFwd}{Mean delta G of the forward primer.}
#'   \item{deltaGRangeFwd}{Range in delta G of the forward primer.}
#'   \item{sequenceFwd}{Sequence of the forward primer, all variants.}
#'   \item{gcContentFwd}{GC-content of the forward primer, all variants.}
#'   \item{tmFwd}{Tm of the forward primer, all variants.}
#'   \item{deltaGFwd}{Delta G of the forward primer, all variants.}
#'   \item{startRev}{Position where the reverse primer starts.}
#'   \item{endRev}{Position where the reverse primer ends.}
#'   \item{lengthRev}{Length of the reverse primer.}
#'   \item{iupacSequenceRev}{IUPAC sequence of the reverse primer.}
#'   \item{coverageRev}{Average coverage of the reverse primer.}
#'   \item{degeneracyRev}{Number of variants of the reverse primer.}
#'   \item{gcContentMeanRev}{Mean GC-content of the reverse primer.}
#'   \item{gcContentRangeRev}{Range in GC-content of the reverse primer.}
#'   \item{tmMeanRev}{Mean tm of the reverse primer.}
#'   \item{tmRangeRev}{Range in tm of the reverse primer.}
#'   \item{deltaGMeanRev}{Mean delta G of the reverse primer.}
#'   \item{deltaGRangeRev}{Range in delta G of the reverse primer.}
#'   \item{sequenceRev}{Sequence of the reverse primer, all variants.}
#'   \item{gcContentRev}{GC-content of the reverse primer, all variants.}
#'   \item{tmRev}{Tm of the reverse primer, all variants.}
#'   \item{deltaGRev}{Delta G of the reverse primer, all variants.}
#'   \item{roiStart}{Start position of the input consensus profile
#'   used for oligo design.}
#'   \item{roiEnd}{End position of the input consensus profile used
#'   for oligo design.}
#' }
#'
#' If a probe is used, the following columns are also included:
#'
#' \describe{
#'   \item{tmDifferencePrimerProbe}{Tm difference between the primer pair
#'   and probe. The tm difference is calculated by subtracting the
#'   mean tm of the probe with the mean tm of the primers.}
#'   \item{deltaGDifferencePrimerProbe}{Difference in delta G between
#'   the primer pair and probe. It is calculated by subtracting the
#'   mean delta G of the probe with the mean delta G of the primers.}
#'   \item{plusPr}{If the probe is valid in positive sense.}
#'   \item{minusPr}{If the probe is valid in negative sense.}
#'   \item{startPr}{Position where the probe starts.}
#'   \item{endPr}{Position where the probe ends.}
#'   \item{lengthPr}{Length of the probe.}
#'   \item{iupacSequencePr}{IUPAC sequence of the probe.}
#'   \item{coveragePr}{Average coverage of the probe.}
#'   \item{degeneracyPr}{Number of variants of the probe.}
#'   \item{gcContentMeanPr}{Mean GC-content of the probe.}
#'   \item{gcContentRangePr}{Range in GC-content of the probe.}
#'   \item{tmMeanPr}{Mean tm of the probe.}
#'   \item{tmRangePr}{Range in tm of the probe.}
#'   \item{deltaGMeanPr}{Mean delta G of the probe.}
#'   \item{deltaGRangePr}{Range in delta G of the probe.}
#'   \item{sequencePr}{Sequence of the probe, all variants.}
#'   \item{gcContentPr}{GC-content of the probe, all variants.}
#'   \item{tmPr}{Tm of the probe, all variants.}
#'   \item{deltaGPr}{Delta G of the probe, all variants.}
#' }
#'
#' @export
#'
#' @examples
#' data("exampleRprimerOligo")
#' ## Get assays using default settings
#' assays(exampleRprimerOligo)
assays <- function(x,
                   lengthRange = c(65, 120),
                   tmDiffPrimers = 5,
                   tmDiffPrimersProbe = c(-10, 10)) {
    if (!methods::is(x, "RprimerOligo")) {
        stop("'x' must be an RprimerOligo object.")
    }
    if (!(min(lengthRange) >= 40 && max(lengthRange) <= 5000)) {
        stop("'lengthRange' must be from 40 to 5000.", call. = FALSE)
    }
    if (!(tmDiffPrimers >= 0 && tmDiffPrimers <= 20)) {
        stop("'tmDiffPrimers' must be from 0 to 20.", call. = FALSE)
    }
    if (!(min(tmDiffPrimersProbe) >= -20 &&
        max(tmDiffPrimersProbe) <= 20)) {
        stop(
            "'tmDiffPrimersProbe' must be from -20 to 20, e.g. c(-1, 5).",
            call. = FALSE
        )
    }
    x <- as.data.frame(x)
    assays <- .combinePrimers(
        x[x$type == "primer", ],
        lengthRange,
        tmDiffPrimers
    )
    if (any(x$type == "probe")) {
        assays <- .addProbes(
            assays,
            x[x$type == "probe", ],
            tmDiffPrimersProbe
        )
    }
    assays <- .beautifyPrimers(assays)
    if (any(x$type == "probe")) {
        assays <- .beautifyProbes(assays)
    }
    RprimerAssay(assays)
}

# Helpers =====================================================================

#' Find all possible primer pairs
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerOligo")
#' x <- exampleRprimerOligo
#' .pairPrimers(x)
.pairPrimers <- function(x) {
    fwd <- x[x$fwd, ]
    rev <- x[x$rev, ]
    pairs <- expand.grid(
        fwd$iupacSequence, rev$iupacSequenceRc, stringsAsFactors = FALSE
    )
    names(pairs) <- c("fwd", "rev")
    fwd <- x[match(pairs$fwd, x$iupacSequence), ]
    rev <- x[match(pairs$rev, x$iupacSequenceRc), ]
    names(fwd) <- paste0(names(fwd), "Fwd")
    names(rev) <- paste0(names(rev), "Rev")
    cbind(fwd, rev)
}

#' Combine primers to assays
#'
#' @inheritParams assays
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerOligo")
#' x <- exampleRprimerOligo
#' .combinePrimers(x)
.combinePrimers <- function(x,
                            lengthRange = c(65, 120),
                            tmDiffPrimers = 2) {
    assays <- .pairPrimers(x)
    ampliconLength <- assays$endRev - assays$startFwd + 1
    tmDifferencePrimer <- abs(assays$tmMeanFwd - assays$tmMeanRev)
    deltaGDifferencePrimer <- abs(
        assays$deltaGMeanFwd - assays$deltaGMeanRev
    )
    start <- assays$startFwd
    end <- assays$endRev
    totalDegeneracy <- assays$degeneracyFwd + assays$degeneracyRev
    score <- assays$scoreFwd + assays$scoreRev
    assays <- cbind(
        start, end, ampliconLength,
        tmDifferencePrimer, deltaGDifferencePrimer, totalDegeneracy, score,
        assays
    )
    assays <- assays[assays$ampliconLength >= min(lengthRange), , drop = FALSE]
    assays <- assays[assays$ampliconLength <= max(lengthRange), , drop = FALSE]
    assays <- assays[assays$tmDifferencePrimer <= tmDiffPrimers, , drop = FALSE]
    if (nrow(assays) == 0L) {
        stop("No assays were found.", call. = FALSE)
    }
    assays
}

#' Find all possible probe candidates to all primer pairs
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerOligo")
#' x <- exampleRprimerOligo
#' assays <- .combinePrimers(x)
#' .identifyProbes(assays, x[x$type == "probe", ])
.identifyProbes <- function(x, probes) {
    lapply(seq_len(nrow(x)), function(i) {
        from <- x$endFwd[[i]] + 1
        to <- x$startRev[[i]] - 1
        probes[probes$start >= from & probes$end <= to, ]
    })
}

#' Combine primers and probes
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerOligo")
#' x <- exampleRprimerOligo
#' probes <- .identifyProbes(assays, x[x$type == "probe", ])
#' test <- .extractProbes(assays, probes)
.extractProbes <- function(x, probes, tmDiffPrimersProbe = c(0, 20)) {
    nProbes <- vapply(probes, nrow, integer(1), USE.NAMES = FALSE)
    select <- lapply(
        seq_along(nProbes), function(x) {
            rep(x, nProbes[[x]])
        }
    )
    select <- unlist(select)
    x <- x[select, , drop = FALSE]
    if (nrow(x) == 0L) {
        stop("No assays with probes could be generated.", call. = FALSE)
    }
    probes <- do.call("rbind", probes)
    names(probes) <- paste0(names(probes), "Pr")
    x <- cbind(x, probes)
    x$totalDegeneracy <- x$totalDegeneracy + probes$degeneracyPr
    x$score <- x$score + probes$scorePr
    tmDifferencePrimerProbe <- vapply(seq_len(nrow(x)), function(i) {
        x$tmMeanPr[[i]] - mean(c(
            x$tmFwd[[i]], x$tmRev[[i]]
        ))
    }, double(1), USE.NAMES = FALSE)
    deltaGDifferencePrimerProbe <- vapply(seq_len(nrow(x)), function(i) {
        x$deltaGMeanPr[[i]] - mean(c(
            x$deltaGFwd[[i]], x$deltaGRev[[i]]
        ))
    }, double(1), USE.NAMES = FALSE)
    x <- cbind(x, tmDifferencePrimerProbe, deltaGDifferencePrimerProbe)
    x <- x[
        x$tmDifferencePrimerProbe >= min(tmDiffPrimersProbe) &
            x$tmDifferencePrimerProbe <= max(tmDiffPrimersProbe),
    ]
    if (nrow(x) == 0L) {
        stop("No assays with probes could be generated.", call. = FALSE)
    }
    x
}

#' Add probes to primer pairs.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerOligo")
#' x <- exampleRprimerOligo
#' assays <- .combinePrimers(x)
#' test <- .addProbes(assays, x[x$type == "probe", ])
.addProbes <- function(x, probes, tmDiffPrimersProbe = c(0, 20)) {
    probeCandidates <- .identifyProbes(x, probes)
    .extractProbes(x, probeCandidates, tmDiffPrimersProbe)
}

.beautifyPrimers <- function(x) {
    drop <- c(
        "typeFwd", "fwdFwd", "revFwd", "typeRev", "fwdRev", "revRev",
        "iupacSequenceRcFwd", "sequenceRcFwd", "iupacSequenceRev",
        "sequenceRev", "scoreFwd", "scoreRev", "roiStartFwd", "roiEndFwd"
    )
    x <- x[!(names(x) %in% drop)]
    names(x)[
        names(x) %in% c("iupacSequenceRcRev", "sequenceRcRev")
    ] <- c("iupacSequenceRev", "sequenceRev")
    names(x)[
        names(x) %in% c("roiStartRev", "roiEndRev")
    ]<- c("roiStart", "roiEnd")
    x <- x[order(x$start), ]
    rownames(x) <- NULL
    x
}

.beautifyProbes <- function(x) {
    drop <- c("typePr", "scorePr", "roiStartPr", "roiEndPr")
    x <- x[!(names(x) %in% drop)]
    rename <- c("fwdPr", "revPr")
    names(x)[names(x) %in% rename] <- c("plusPr", "minusPr")
    moveToLast <- c("roiStart", "roiEnd")
    x <- x[c(setdiff(names(x), moveToLast), moveToLast)]
    moveToFirst <- c(
        "start", "end", "ampliconLength", "tmDifferencePrimer",
        "tmDifferencePrimerProbe"
    )
    x <- x[c(moveToFirst, setdiff(names(x), moveToFirst))]
    x
}
