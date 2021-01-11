#' Design (RT)-PCR assays
#'
#' \code{assays()} combines forward and reverse primers
#' and probes (if selected) to (RT)-PCR assays.
#'
#' @param x An \code{RprimerOligo} object, with or without probes.
#'
#' @param lengthRange
#' Range of the amplicon range, a numeric vector [40, 5000], defaults to
#' \code{c(65, 120)}.
#'
#' @param tmDiffPrimers
#' Maximum Tm difference between the two primers
#' (absolute value). A number [0, 30], defaults to 2.
#'
#' @param tmDiffPrimersProbe
#' Acceptable Tm difference between the primers (average Tm of the
#' primer pair) and probe. A numeric vector [-20, 20],
#' defaults to \code{c(0, 20)}.
#' The Tm-difference is calculated by subtracting the
#' Tm of the probe with the average Tm of the majority
#' primer pair. Thus, a negative Tm-difference
#' means that the Tm of the probe is lower than the average Tm of the
#' primer pair.
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
#'   \item{tmDifferencePrimer}{Difference in Tm between
#'   the forward and reverse primer, absolute value.}
#'   \item{totalDegeneracy}{Total number of oligos in the assay.}
#'   \item{startFwd}{Position where the forward primer starts.}
#'   \item{endFwd}{Position where the reverse primer ends.}
#'   \item{lengthFwd}{Length of the forward primer.}
#'   \item{majorityFwd}{Majority sequence of the forward primer.}
#'   \item{gcMajorityFwd}{GC-content of the forward primer
#'   (majority sequence), proportion.}
#'   \item{identityFwd}{Average identity of the forward primer.}
#'   \item{tmMajorityFwd}{Tm of the forward primer
#'   (majority sequence).}
#'   \item{iupacFwd}{IUPAC sequence (i.e. with degenerate bases)
#'   of the forward primer.}
#'   \item{degeneracyFwd}{Number of variants of the forward primer}.
#'   \item{startRev}{Position where the reverse primer starts.}
#'   \item{endRev}{Position where the reverse primer ends.}
#'   \item{lengthRev}{Length of the reverse primer.}
#'   \item{majorityRev}{Majority sequence of the reverse primer.}
#'   \item{gcMajorityRev}{GC-content of the reverse primer
#'   (majority sequence), proportion.}
#'   \item{tmMajorityRev}{Tm of the reverse primer
#'   (majority sequence).}
#'   \item{identityRev}{Average identity of the reverse primer.}
#'   \item{iupacRev}{IUPAC sequence (i.e. with degenerate bases)
#'   of the reverse primer.}
#'   \item{degeneracyRev}{Number of variants of the reverse primer.}
#'   \item{allFwd}{Lists with all sequence variants of the forward primer.}
#'   \item{gcAllFwd}{Lists with the GC content of all
#'   sequence variants of the forward primer.}
#'   \item{tmAllFwd}{Lists with the Tm of all sequence variants of
#'   the reverse primer.}
#'   \item{allRev}{Lists with all sequence variants of the reverse primer.}
#'   \item{gcAllRev}{Lists with the GC content of all
#'   sequence variants of the forward primer.}
#'   \item{tmAllRev}{Lists with the Tm of all sequence variants of
#'   the reverse primer.}
#'   \item{alignmentStart}{Start position of the input consensus profile
#'   used for oligo design.}
#'   \item{alingnmentEnd}{End position of the input consensus profile used
#'   for oligo design.}
#' }
#'
#' If a probe is used, the following columns are also included:
#'
#' \describe{
#'   \item{tmDifferencePrimerProbe}{Difference in Tm between the average
#'   Tm of the primer pair and the probe, majority sequences.}
#'   \item{startPr}{Position where the probe starts.}
#'   \item{endPr}{Position where the probe ends.}
#'   \item{lengthPr}{Length of the probe.}
#'   \item{majorityPr}{Majority sequence of the probe.}
#'   \item{gcMajorityPr}{GC-content of the probe
#'   (majority sequence), proportion.}
#'   \item{tmMajorityPr}{Tm of the probe
#'   (majority sequence).}
#'   \item{identityPr}{Average identity of the probe.}
#'   \item{iupacPr}{IUPAC sequence (i.e. with degenerate bases)
#'   of the probe.}
#'   \item{degeneracyPr}{Number of variants of the probe.}
#'   \item{sensePr}{Sense of the probe (pos or neg). If both probes are valid,
#'   the probe with the least G:s is selected.}
#'   \item{allPr}{Lists with all sequence variants of the probe.}
#'   \item{gcAllPr}{Lists with the GC content of all
#'   sequence variants of the probe.}
#'   \item{tmAllPr}{Lists with the Tm of all sequence variants of
#'   the probe.}
#' }
#'
#' @seealso getOligos
#'
#' @examples
#' data("exampleRprimerOligo")
#' ## Get assays using default settings
#' assays(exampleRprimerOligo)
#' @export
assays <- function(x,
                   lengthRange = c(65, 120),
                   tmDiffPrimers = 2,
                   tmDiffPrimersProbe = c(0, 20)) {
    if (!methods::is(x, "RprimerOligo")) {
        stop("'x' must be an RprimerOligo object.")
    }
    if (!(min(lengthRange) >= 40 && max(lengthRange) <= 5000)) {
        stop("'lengthRange' must be from 40 to 5000.", call. = FALSE)
    }
    if (!(tmDiffPrimers >= 0 && tmDiffPrimers <= 30)) {
        stop("'tmDiffPrimers' must be from 0 to 30.", call. = FALSE)
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
    # assays <- .beautifyAssays(assay.s)
    RprimerAssay(assays)
}

# Helpers =====================================================================

#' Find all possible primer pairs
#'
#' @keywords internal
#'
#' @noRd
.pairPrimers <- function(x) {
    x <- as.data.frame(x)
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

#' Combine primers to assays
#'
#' @inheritParams assays
#'
#' @keywords internal
#'
#' @noRd
.combinePrimers <- function(x,
                            lengthRange = c(65, 120),
                            tmDiffPrimers = 2) {
    assays <- .pairPrimers(x)
    ampliconLength <- assays$endRev - assays$startFwd + 1
    tmDifferencePrimer <- abs(assays$tmMeanFwd - assays$tmMeanRev)
    start <- assays$startFwd
    end <- assays$endRev
    totalDegeneracy <- assays$degeneracyFwd + assays$degeneracyRev
    assays <- cbind(
        start, end, ampliconLength,
        tmDifferencePrimer, totalDegeneracy,
        assays
    )
    assays <- assays[assays$ampliconLength >= min(lengthRange), ]
    assays <- assays[assays$ampliconLength <= max(lengthRange), ]
    assays <- assays[assays$tmDifferencePrimer <= tmDiffPrimers, ]
    if (nrow(assays) == 0) {
        stop("No assays were found.", call. = FALSE)
    }
    assays
}

#' Find all possible probe candidates to all primer pairs
#'
#' @keywords internal
#'
#' @noRd
.identifyProbes <- function(x, probes) {
    lapply(seq_len(nrow(x)), function(i) {
        from <- x$endFwd[[i]] + 1
        to <- x$startRev[[i]] - 1
        probes[probes$start >= from & probes$end <= to, ]
    })
}

#' Combine primers and probes
#'
#' @inheritParams assays
#'
#' @keywords internal
#'
#' @noRd
.extractProbes <- function(x, probes, tmDiffPrimersProbe = 2) {
    nProbes <- vapply(probes, nrow, integer(1), USE.NAMES = FALSE)
    select <- lapply(
        seq_along(nProbes), function(x) {
            rep(x, nProbes[[x]])
        }
    )
    select <- unlist(select)
    x <- x[select, , drop = FALSE]
    if (nrow(x) == 0) {
        stop("No assays with probes could be generated.", call. = FALSE)
    }
    probes <- do.call("rbind", probes)
    names(probes) <- paste0(names(probes), "Pr")
    x <- cbind(x, probes)
    x$totalDegeneracy <- x$totalDegeneracy + probes$degeneracyPr
    tmDifferencePrimerProbe <- vapply(seq_len(nrow(x)), function(i) {
        x$tmMeanPr[[i]] - mean(c(
            x$tmMeanFwd[[i]], x$tmMeanRev[[i]]
        ))
    }, double(1), USE.NAMES = FALSE)
    x <- cbind(x, tmDifferencePrimerProbe)
    x <- x[
        x$tmDifferencePrimerProbe >= min(tmDiffPrimersProbe) &
            x$tmDifferencePrimerProbe <= max(tmDiffPrimersProbe),
    ]
    x
}

#' Add probes to primer pairs.
#'
#' @param x Assays to add probes to.
#'
#' @param probes Candidate probes.
#'
#' @inheritParams assays
#'
#' @return A data frame.
#'
#' @keywords internal
#'
#' @noRd
.addProbes <- function(x, probes, tmDiffPrimersProbe = c(0, 20)) {
    probeCandidates <- .identifyProbes(x, probes)
    assays <- .extractProbes(x, probeCandidates, tmDiffPrimersProbe)
    if (nrow(assays) == 0) {
        stop("No assays with probes could be generated.", call. = FALSE)
    }
    assays
}

.beautifyAssays <- function(x) {
    drop <- c(
        "fwdFwd", "revFwd", "fwdRev", "revRev", "typeFwd", "typeRev",
        "iupacSequenceRcFwd", "sequenceRcFwd", "iupacSequenceRev",
        "sequenceRev"
    )
    assays <- assays[!(names(assays) %in% drop)]
    names(assays)[grep("Rc", names(assays))] <- gsub(
        "Rc", "", names(assays)[grep("Rc", names(assays))]
    )
    # beautify
    #     drop <- c("majorityRc", "iupacRc")
    # probes <- probes[!names(probes) %in% drop]
    rownames(x) <- NULL
    alignmentStart <- assays$alignmentStartFwd
    alignmentEnd <- assays$alignmentEndFwd
    assays <- assays[-grep("type", names(assays))]
    assays <- assays[-grep("alignmentStart", names(assays))]
    assays <- assays[-grep("alignmentEnd", names(assays))]
    assays <- cbind(assays, alignmentStart, alignmentEnd)
    assays <- assays[order(assays$start), ]
}
