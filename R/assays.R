#' Design (RT)-PCR assays
#'
#' \code{assays()} combines forward primers, reverse primers
#' and probes to (RT)-PCR assays from an \code{RprimerOligo} object.
#'
#' @param x
#' An \code{RprimerOligo} object, which can be with or without probes.
#'
#' @param length
#' Amplicon length range, a numeric vector [40, 5000], defaults to
#' \code{c(65, 120)}.
#'
#' @param tmDifferencePrimers
#' Maximum allowed difference between the mean tm of the forward and reverse
#' primer (in Celcius degrees, as an absolute value). Defaults to \code{NULL},
#' which means that primers will be paired regardless of their tm.
#'
#' @return
#' An \code{RprimerAssay} object, which contains the following information:
#'
#' \describe{
#'   \item{start}{Position where the assay starts.}
#'   \item{end}{Position where the assay ends.}
#'   \item{length}{Length of the amplicon.}
#'   \item{totalDegeneracy}{Total number of oligos in the assay.}
#'   \item{score}{Summarized oligo score. The lowest, and best,
#'   possible score is 0. The highest possible score is 24 for assays
#'   with only primers,
#'   and 36 for assays with probes.
#'   See \code{?oligos} for more information about the scoring system.}
#'   \item{startFwd}{Start position of the forward primer.}
#'   \item{endFwd}{End positon of the forward primer.}
#'   \item{lengthFwd}{Length of the forward primer.}
#'   \item{iupacSequenceFwd}{Forward primer sequence in IUPAC format
#'   (i.e. with ambiguous bases).}
#'   \item{identityFwd}{For ambiguous primers: Average identity of the
#'   forward primer. For mixed primers: Average identity of the 5' (consensus)
#'   part of the forward primer. The value can range from 0 to 1.}
#'   \item{coverageFwd}{For ambiguous primers: Average coverage of the
#'   forward primer. For mixed primers: Average coverage of the 3' (degenerate)
#'   part of the forward primer. The value can range from 0 to 1.}
#'   \item{degeneracyFwd}{Number of sequence variants of the forward primer.}
#'   \item{gcContentMeanFwd}{Mean GC-content of all sequence variants of the
#'   forward primer.}
#'   \item{gcContentRangeFwd}{Range in GC-content of all sequence variants of
#'     the forward primer.}
#'   \item{tmMeanFwd}{Mean tm of all sequence variants of the forward primer
#'   (in Celcius degrees).}
#'   \item{tmRangeFwd}{Range in tm of all sequence variants of the forward
#'   primer (in Celcius degrees).}
#'   \item{deltaGMeanFwd}{Mean delta G of all sequence variants of the
#'   forward primer (in kcal/mol).}
#'   \item{deltaGRangeFwd}{Range in delta G of all sequence variants of the
#'   forward primer (in kcal/mol).}
#'   \item{sequenceFwd}{All sequence variants of the forward primer.}
#'   \item{gcContentFwd}{GC-content of all sequence variants of the forward
#'   primer.}
#'   \item{tmFwd}{Tm of all sequence variants of the forward primer
#'   (in Celcius degrees).}
#'   \item{deltaGFwd}{Delta G of all sequence variants
#'   of the forward primer (in kcal/mol).}
#'   \item{methodFwd}{Design method used to generate the forward
#'   primer: "ambiguous" or "mixedFwd.}
#'   \item{startRev}{Start position of the reverse primer.}
#'   \item{endRev}{End positon of the reverse primer.}
#'   \item{lengthRev}{Length of the reverse primer.}
#'   \item{iupacSequenceRev}{Reverse primer sequence in IUPAC format
#'   (i.e. with ambiguous bases).}
#'   \item{identityRev}{For ambiguous primers: Average identity of the
#'   reverse primer. For mixed primers: Average identity of the 5' (consensus)
#'   part of the reverse primer. The value can range from 0 to 1.}
#'   \item{coverageRev}{For ambiguous primers: Average coverage of the
#'   reverse primer. For mixed primers: Average coverage of the 3' (degenerate)
#'   part of the reverse primer. The value can range from 0 to 1.}
#'   \item{degeneracyRev}{Number of sequence variants of the reverse primer.}
#'   \item{gcContentMeanRev}{Mean GC-content of all sequence variants of the
#'   reverse primer.}
#'   \item{gcContentRangeRev}{Range in GC-content of all sequence variants of
#'     the reverse primer.}
#'   \item{tmMeanRev}{Mean tm of all sequence variants of the reverse primer
#'   (in Celcius degrees).}
#'   \item{tmRangeRev}{Range in tm of all sequence variants of the reverse
#'   primer (in Celcius degrees).}
#'   \item{deltaGMeanRev}{Mean delta G of all sequence variants of the
#'   reverse primer (in kcal/mol).}
#'   \item{deltaGRangeRev}{Range in delta G of all sequence variants of the
#'   reverse primer (in kcal/mol).}
#'   \item{sequenceRev}{All sequence variants of the reverse primer.}
#'   \item{gcContentRev}{GC-content of all sequence variants of the reverse
#'   primer.}
#'   \item{tmRev}{Tm of all sequence variants of the reverse primer
#'   (in Celcius degrees).}
#'   \item{deltaGRev}{Delta G of all sequence variants
#'   of the reverse primer (in kcal/mol).}
#'   \item{methodRev}{Design method used to generate the forward
#'   primer: "ambiguous" or "mixedRev.}
#'   \item{roiStart}{Start position of the input consensus profile
#'   used for oligo design.}
#'   \item{roiEnd}{End position of the input consensus profile used
#'   for oligo design.}
#' }
#'
#' If a probe is included in the input \code{RprimerOligo} object,
#' the following columns are also included:
#'
#' \describe{
#'   \item{startPr}{Start position of the probe.}
#'   \item{endPr}{End positon of the probe.}
#'   \item{lengthPr}{Length of the probe.}
#'   \item{iupacSequencePr}{Probe sequence in plus sense, in IUPAC format.}
#'   \item{iupacSequenceRcPr}{Probe sequence in minus sense,
#'   in IUPAC format.}
#'   \item{identityPr}{For ambiguous primers: Average identity of the
#'   probe. For mixed primers: Average identity of the 5' (consensus)
#'   part of the probe. The value can range from 0 to 1.}
#'   \item{coveragePr}{For ambiguous primers: Average coverage of the
#'   probe. For mixed primers: Average coverage of the 3' (degenerate)
#'   part of the probe. The value can range from 0 to 1.}
#'   \item{degeneracyPr}{Number of sequence variants of the probe.}
#'   \item{gcContentMeanPr}{Mean GC-content of all sequence variants of the
#'   probe.}
#'   \item{gcContentRangePr}{Range in GC-content of all sequence variants of
#'     the probe.}
#'   \item{tmMeanPr}{Mean tm of all sequence variants of the probe
#'   (in Celcius degrees).}
#'   \item{tmRangePr}{Range in tm of all sequence variants of the forward
#'   primer (in Celcius degrees).}
#'   \item{deltaGMeanPr}{Mean delta G of all sequence variants of the
#'   probe (in kcal/mol).}
#'   \item{deltaGRangePr}{Range in delta G of all sequence variants of the
#'   probe (in kcal/mol).}
#'   \item{sequencePr}{All sequence variants of the probe, in plus sense.}
#'   \item{sequenceRcPr}{All sequence variants of the probe, in minus sense.}
#'   \item{gcContentPr}{GC-content of all sequence variants of the probe.}
#'   \item{tmPr}{Tm of all sequence variants of the probe
#'   (in Celcius degrees).}
#'   \item{deltaGPr}{Delta G of all sequence variants
#'   of the probe (in kcal/mol).}
#'   \item{methodPr}{Design method used to generate the probe.}
#'   \item{plusPr}{If the probe is valid in plus sense.}
#'   \item{minusPr}{If the probe is valid in minus sense.}
#' }
#'
#' An error message will return if no assays are found.
#'
#' @export
#'
#' @examples
#' data("exampleRprimerOligo")
#'
#' ## Design assays using default settings
#' assays(exampleRprimerOligo)
#'
#' ## Modify the length range
#' assays(exampleRprimerOligo, length = c(1000, 2000))
assays <- function(x,
                   length = c(65, 120),
                   tmDifferencePrimers = NULL) {
    if (!methods::is(x, "RprimerOligo")) {
        stop("'x' must be an RprimerOligo object.")
    }
    if (!(min(length) >= 40 && max(length) <= 5000)) {
        stop("'length' must be from 40 to 5000.", call. = FALSE)
    }
    if (!is.null(tmDifferencePrimers) && !is.numeric(tmDifferencePrimers)) {
        stop("'tmDifferencePrimers must be either 'NULL' or a number.")
    }
    x <- as.data.frame(x)
    assays <- .combinePrimers(x[x$type == "primer", ], length)
    if (any(x$type == "probe")) {
        assays <- .addProbes(assays, x[x$type == "probe", ])
    }
    assays <- .beautifyPrimers(assays)
    if (any(x$type == "probe")) {
        assays <- .beautifyProbes(assays)
    }
    RprimerAssay(assays)
}

# Helpers ======================================================================

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
        fwd$iupacSequence, rev$iupacSequenceRc,
        stringsAsFactors = FALSE
    )
    names(pairs) <- c("fwd", "rev")
    fwd <- x[match(pairs$fwd, x$iupacSequence), ]
    rev <- x[match(pairs$rev, x$iupacSequenceRc), ]
    names(fwd) <- paste0(names(fwd), "Fwd")
    names(rev) <- paste0(names(rev), "Rev")
    cbind(fwd, rev)
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerOligo")
#' x <- exampleRprimerOligo
#' .combinePrimers(x)
.combinePrimers <- function(x,
                            length = c(65, 120),
                            tmDifferencePrimers = NULL) {
    assays <- .pairPrimers(x)
    ampLength <- assays$endRev - assays$startFwd + 1
    start <- assays$startFwd
    end <- assays$endRev
    totalDegeneracy <- assays$degeneracyFwd + assays$degeneracyRev
    score <- assays$scoreFwd + assays$scoreRev
    assays <- cbind(
        start, end, "length" = ampLength, totalDegeneracy, score,
        assays
    )
    assays <- assays[assays$length >= min(length), , drop = FALSE]
    assays <- assays[assays$length <= max(length), , drop = FALSE]
    if (!is.null(tmDifferencePrimers)) {
        tmDifferencePrimers <- abs(tmDifferencePrimers)
        tmDifference <- abs(assays$tmMeanFwd - assays$tmMeanRev)
        assays <- assays[tmDifference <= tmDifferencePrimers, , drop = FALSE]
    }
    if (nrow(assays) == 0L) {
        stop("No assays were found.", call. = FALSE)
    }
    assays
}

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

#' @noRd
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
    ] <- c("roiStart", "roiEnd")
    x <- x[order(x$start), ]
    rownames(x) <- NULL
    x
}

#' @noRd
.beautifyProbes <- function(x) {
    drop <- c("typePr", "scorePr", "roiStartPr", "roiEndPr")
    x <- x[!(names(x) %in% drop)]
    rename <- c("fwdPr", "revPr")
    names(x)[names(x) %in% rename] <- c("plusPr", "minusPr")
    moveToLast <- c("roiStart", "roiEnd")
    x <- x[c(setdiff(names(x), moveToLast), moveToLast)]
    moveToFirst <- c("start", "end", "length")
    x <- x[c(moveToFirst, setdiff(names(x), moveToFirst))]
    x
}
