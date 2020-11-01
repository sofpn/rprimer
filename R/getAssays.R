#' Get (RT)-PCR assays from oligos
#'
#' \code{getAssays()} combines forward and reverse primers
#' and (if selected) probes to (RT)-PCR assays.
#'
#' @param x An \code{RprimerOligo} object, with or without probes.
#'
#' @param length
#' Amplicon length, a numeric vector [40, 5000]. Defaults to
#' \code{65:120}.
#'
#' @param maxTmDifferencePrimers
#' Maximum Tm difference between the two primers
#' (absolute value). A number [0, 30]. Defaults to 2.
#' Note that the Tm-difference is calculated from the majority primers, and
#' may thus be misleading for degenerate (IUPAC) primers.
#'
#' @param tmDifferencePrimersProbe
#' Acceptable Tm difference between the primers (average Tm of the
#' primer pair) and probe. A numeric vector [-20, 20],
#' defaults to \code{c(0, 20)}.
#' The Tm-difference is calculated by subtracting the
#' Tm of the probe with the average Tm of the majority
#' primer pair. Thus, a negative Tm-difference
#' means that the Tm of the probe is lower than the average Tm of the
#' primer pair.
#' Note that the Tm-difference is calculated from the majority oligos, and
#' may thus be misleading for degenerate (IUPAC) oligos.
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
#'   \item{meanIdentity}{Average identity score of the primers
#'   (and probe if selected)}.
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
#' }
#'
#' If the option \code{showAllVariants == TRUE} was used for primer design
#' in \code{get_oligos()}, the following columns are also added:
#'
#' \describe{
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
#' }
#'
#' And, if the option \code{showAllVariants == TRUE} was used
#' in \code{getOligos()}, the following data are also added:
#'
#' \describe{
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
#' getAssays(exampleRprimerOligo)
#'
#' @export
getAssays <- function(x,
                      length = 65:120,
                      maxTmDifferencePrimers = 2,
                      tmDifferencePrimersProbe = c(0, 20)) {
    if (!methods::is(x, "RprimerOligo")) {
        stop("'x' must be an RprimerOligo object.")
    }
    x <- as.data.frame(x)
    assays <- .combinePrimers(
        primers = x[x$type == "primer", ], length = length,
        maxTmDifferencePrimers = maxTmDifferencePrimers
    )
    if (any(x$type == "probe")) {
        assays <- .addProbes(
            assays = assays, probes = x[x$type == "probe", ],
            tmDifferencePrimersProbe = tmDifferencePrimersProbe
        )
    }
    alignmentStart <- assays$alignmentStartFwd
    alignmentEnd <- assays$alignmentEndFwd
    assays <- assays[-grep("type", names(assays))]
    assays <- assays[-grep("alignmentStart", names(assays))]
    assays <- assays[-grep("alignmentEnd", names(assays))]
    assays <- tibble::add_column(assays, alignmentStart, alignmentEnd)
    assays <- assays[order(assays$start), ]
    RprimerAssay(assays)
}

# Helpers =====================================================================

#' Calculate G content of a DNA sequence
#'
#' \code{.gContent()} finds the G content of a DNA sequence.
#'
#' @param x
#' A DNA sequence (a character vector of length one).
#'
#' @return The G content of x.
#'
#' @keywords internal
#'
#' @noRd
.gContent <- function(x) {
    x <- .splitSequence(x)
    gCount <- length(which(x == "G"))
    totalCount <- length(which(x == "A" | x == "C" | x == "G" | x == "T"))
    gCount / totalCount
}

#' Combine primers to assays
#'
#' @inheritParams getAssays
#'
#' @keywords internal
#'
#' @noRd
.combinePrimers <- function(primers,
                            length = 65:120,
                            maxTmDifferencePrimers = 2) {
    if (!(maxTmDifferencePrimers > 0 && maxTmDifferencePrimers < 30)) {
        stop("'maxTmDifferencePrimers' must be from 0 to 30.", call. = FALSE)
    }
    if (!(min(length) >= 40 && max(length) <= 5000)) {
        stop("'length' must be from 40 to 5000.", call. = FALSE)
    }
    fwd <- primers[!is.na(primers$majority), ]
    rev <- primers[!is.na(primers$majorityRc), ]
    combinations <- expand.grid(
        fwd$majority, rev$majorityRc,
        stringsAsFactors = FALSE
    )
    names(combinations) <- c("fwdMajority", "revMajority")
    indexFwd <- match(combinations$fwdMajority, primers$majority)
    indexRev <- match(combinations$revMajority, primers$majorityRc)
    fwd <- primers[indexFwd, ]
    rev <- primers[indexRev, ]
    names(fwd) <- paste0(names(fwd), "Fwd")
    names(rev) <- paste0(names(rev), "Rev")
    assays <- dplyr::bind_cols(fwd, rev)
    assays <- tibble::tibble(assays)
    ampliconLength <- assays$endRev - assays$startFwd + 1
    ampliconLength <- as.integer(ampliconLength)
    tmDifferencePrimer <- abs(assays$tmMajorityFwd - assays$tmMajorityRev)
    start <- assays$startFwd
    end <- assays$endRev
    meanIdentity <- mean(c(assays$identityFwd, assays$identityRev))
    totalDegeneracy <- assays$degeneracyFwd + assays$degeneracyRev
    assays <- tibble::add_column(
        assays, start, end, ampliconLength,
        tmDifferencePrimer, meanIdentity, totalDegeneracy,
        .before = "startFwd"
    )
    drop <- c("majorityRcFwd", "iupacRcFwd", "majorityRev", "iupacRev")
    assays <- assays[!(names(assays) %in% drop)]
    if (any(grep("all", names(assays)))) {
        drop <- c("allRcFwd", "allRev")
        assays <- assays[!(names(assays) %in% drop)]
    }
    names(assays)[grep("Rc", names(assays))] <- gsub(
        "Rc", "", names(assays)[grep("Rc", names(assays))]
    )
    assays <- assays[assays$ampliconLength >= min(length), ]
    assays <- assays[assays$ampliconLength <= max(length), ]
    assays <- assays[assays$tmDifferencePrimer <= maxTmDifferencePrimers, ]
    if (nrow(assays) == 0L) {
        stop("No assays were found.", call. = FALSE)
    }
    assays
}

#' Add probes to (RT)-PCR assays
#'
#' \code{.addProbes()} adds probes to (RT)-PCR assays.
#'
#' @param assays Assays to add probes to.
#'
#' @param probes Candidate probes.
#'
#' @inheritParams getAssays
#'
#' @return Assays with probes. A tibble (a data frame).
#'
#' @keywords internal
#'
#' @noRd
.addProbes <- function(assays, probes, tmDifferencePrimersProbe = c(0, 20)) {
    if (!(is.numeric(tmDifferencePrimersProbe) &&
        min(tmDifferencePrimersProbe) >= -20 &&
        max(tmDifferencePrimersProbe) <= 20)
    ) {
        stop(
            "'tmDifference' must be from -20 to 20, e.g. c(-1, 5).",
            call. = FALSE
        )
    }
    probeCandidates <- purrr::map(seq_len(nrow(assays)), function(i) {
        from <- assays$endFwd[[i]] + 2
        to <- assays$startRev[[i]] - 2
        probe <- probes[probes$start >= from & probes$end <= to, ]
        probe
    })
    numberOfProbes <- purrr::map_int(probeCandidates, nrow)
    rowsToSelect <- purrr::map(
        seq_along(numberOfProbes), ~ rep(.x, numberOfProbes[[.x]])
    )
    rowsToSelect <- unlist(rowsToSelect, use.names = FALSE)
    assays <- assays[rowsToSelect, ]
    if (nrow(assays) == 0L) {
        stop("No assays with probes could be generated.", call. = FALSE)
    }
    probeCandidates <- do.call("rbind", probeCandidates)
    sense <- purrr::map2_chr(
        probeCandidates$majority, probeCandidates$majorityRc, function(x, y) {
            if (is.na(x)) sense <- "neg"
            if (is.na(y)) sense <- "pos"
            if (!is.na(x) && !is.na(y)) {
                gContentPos <- .gContent(x)
                gContentNeg <- .gContent(y)
                sense <- ifelse(gContentPos <= gContentNeg, "pos", "neg")
            }
            sense
        }
    )
    probeCandidates$majority <- ifelse(
        sense == "pos", probeCandidates$majority, probeCandidates$majorityRc
    )
    probeCandidates$iupac <- ifelse(
        sense == "pos", probeCandidates$iupac, probeCandidates$iupacRc
    )
    probeCandidates <- tibble::add_column(probeCandidates, sense)
    drop <- c("majorityRc", "iupacRc")
    probeCandidates <- probeCandidates[!names(probeCandidates) %in% drop]
    names(probeCandidates) <- paste0(names(probeCandidates), "Pr")
    assays <- dplyr::bind_cols(assays, probeCandidates)
    assays$meanIdentity <- mean(
        c(assays$identityFwd, assays$identityRev, assays$identityPr)
    )
    assays$totalDegeneracy <- assays$totalDegeneracy + probeCandidates$degeneracyPr
    tmDifferencePrimerProbe <- purrr::map_dbl(
        seq_len(nrow(assays)), function(x) {
            assays$tmMajorityPr[[x]] - mean(
                assays$tmMajorityFwd[[x]], assays$tmMajorityRev[[x]]
            )
        }
    )
    assays <- tibble::add_column(
        assays, tmDifferencePrimerProbe, .after = "tmDifferencePrimer"
    )
    assays <- assays[
        assays$tmDifferencePrimerProbe >= min(tmDifferencePrimersProbe),
    ]
    assays <- assays[
        assays$tmDifferencePrimerProbe <= max(tmDifferencePrimersProbe),
    ]
    if (nrow(assays) == 0L) {
        stop("No assays with probes could be generated.", call. = FALSE)
    }
    assays
}
