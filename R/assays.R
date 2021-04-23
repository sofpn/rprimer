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
#'   \item{sequenceFwd}{Sequence of the forward primer, all variants.}
#'   \item{gcContentFwd}{GC-content of the forward primer, all variants.}
#'   \item{tmFwd}{Tm of the forward primer, all variants.}
#'   \item{dHFwd}{delta H of all sequence variants of the forward
#'   primer (in cal/mol).}
#'   \item{dSFwd}{delta S of all sequence
#'   variants of the forward primer (in cal/K/mol).}
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
#'   \item{sequenceRev}{Sequence of the reverse primer, all variants.}
#'   \item{gcContentRev}{GC-content of the reverse primer, all variants.}
#'   \item{tmRev}{Tm of the reverse primer, all variants.}
#'   \item{dHRev}{delta H of all sequence variants of the reverse
#'   primer (in cal/mol).}
#'   \item{dSRev}{delta S of all sequence
#'   variants of the reverse primer (in cal/K/mol).}
#'   \item{roiStart}{Start position of the input consensus profile
#'   used for oligo design.}
#'   \item{roiEnd}{End position of the input consensus profile used
#'   for oligo design.}
#' }
#'
#' If a probe is used, the following columns are also included:
#'
#' \describe{
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
#'   \item{sequencePr}{Sequence of the probe, all variants.}
#'   \item{gcContentPr}{GC-content of the probe, all variants.}
#'   \item{tmPr}{Tm of the probe, all variants.}
#'   \item{dHPr}{delta H of all sequence variants of the forward
#'   primer (in cal/mol).}
#'   \item{dSPr}{delta S of all sequence
#'   variants of the forward primer (in cal/K/mol).}
#' }
#'
#' @export
#'
#' @examples
#' data("exampleRprimerOligo")
#'
#' ## Design assays using default settings
#' assays(exampleRprimerOligo)
assays <- function(x, lengthRange = c(65, 120)) {
    if (!methods::is(x, "RprimerOligo")) {
        stop("'x' must be an RprimerOligo object.")
    }
    if (!(min(lengthRange) >= 40 && max(lengthRange) <= 5000)) {
        stop("'lengthRange' must be from 40 to 5000.", call. = FALSE)
    }
    x <- as.data.frame(x)
    assays <- .combinePrimers(x[x$type == "primer", ], lengthRange)
    if (any(x$type == "probe")) {
        assays <- .addProbes(assays, x[x$type == "probe", ])
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
.combinePrimers <- function(x, lengthRange = c(65, 120)) {
    assays <- .pairPrimers(x)
    ampliconLength <- assays$endRev - assays$startFwd + 1
    start <- assays$startFwd
    end <- assays$endRev
    totalDegeneracy <- assays$degeneracyFwd + assays$degeneracyRev
    score <- assays$scoreFwd + assays$scoreRev
    assays <- cbind(
        start, end, ampliconLength, totalDegeneracy, score,
        assays
    )
    assays <- assays[assays$ampliconLength >= min(lengthRange), , drop = FALSE]
    assays <- assays[assays$ampliconLength <= max(lengthRange), , drop = FALSE]
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
.extractProbes <- function(x, probes) {
    nProbes <- vapply(probes, nrow, integer(1), USE.NAMES = FALSE)
    select <- lapply(seq_along(nProbes), function(x) {
            rep(x, nProbes[[x]])
    })
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
.addProbes <- function(x, probes) {
    probeCandidates <- .identifyProbes(x, probes)
    .extractProbes(x, probeCandidates)
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
    moveToFirst <- c("start", "end", "ampliconLength")
    x <- x[c(moveToFirst, setdiff(names(x), moveToFirst))]
    x
}
