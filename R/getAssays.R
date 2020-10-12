## exempel ## tester ## s4 class

#' Get (RT)-PCR assays from oligos
#'
#' \code{getAssays()} combines forward and reverse primers
#' and (if selected) probes to (RT)-PCR assays.
#'
#' @param primers An \code{RprimerOligo} object.
#'
#' @param probes
#' Optional, defaults to \code{NULL}.
#' An \code{RprimerOligo} object if probes are to be used.
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
#' @param tmDifferenceProbes
#' Optional. Acceptable Tm difference between the primers (average Tm of the
#' primer pair) and probe. A numeric vector [-20, 20],
#' defaults to \code{c(0, 20)}.
#' The Tm-difference is calculated by subtracting the
#' Tm of the probe with the average Tm of the
#' primer pair. Thus, a negative Tm-difference
#' means that the Tm of the probe is lower than the average Tm of the
#' primer pair.
#' Note that the Tm-difference is calculated from the majority oligos, and
#' may thus be misleading for degenerate (IUPAC) oligos.
#'
#' @return
#' A tibble (a data frame) with all candidate
#' assays. An error message will return if no assays are found.
#'
#' @note Results:
#'
#' The tibble contains the following information:
#'
#' \describe{
#'   \item{Begin}{Position where the assay begins.}
#'   \item{End}{Position where the assay ends.}
#'   \item{Amplicon_length}{Length of the amplicon.}
#'   \item{Tm_difference_primer}{Difference in Tm between
#'   the forward and reverse primer, absolute value.}
#'   \item{Mean_identity}{Average identity score of the primers
#'   (and probe if selected)}.
#'   \item{Total_degeneracy}{Total number of oligos in the assay.}
#'   \item{Begin_fwd}{Position where the forward primer begins.}
#'   \item{End_fwd}{Position where the reverse primer ends.}
#'   \item{Length_fwd}{Length of the forward primer.}
#'   \item{Majority_fwd}{Majority sequence of the forward primer.}
#'   \item{GC_majority_fwd}{GC-content of the forward primer
#'   (majority sequence), proportion.}
#'   \item{Identity_fwd}{Average identity of the forward primer.}
#'   \item{Tm_majority_fwd}{Tm of the forward primer
#'   (majority sequence).}
#'   \item{IUPAC_fwd}{IUPAC sequence (i.e. with degenerate bases)
#'   of the forward primer.}
#'   \item{Degeneracy_fwd}{Number of variants of the forward primer.
#'   \item{Begin_rev}{Position where the reverse primer begins.}
#'   \item{End_rev}{Position where the reverse primer ends.}
#'   \item{Length_rev}{Length of the reverse primer.}
#'   \item{Majority_rev}{Majority sequence of the reverse primer.}
#'   \item{GC_majority_rev}{GC-content of the reverse primer
#'   (majority sequence), proportion.}
#'   \item{Tm_majority_rev}{Tm of the reverse primer
#'   (majority sequence).}
#'   \item{Identity_rev}{Average identity of the reverse primer.}
#'   \item{IUPAC_rev}{IUPAC sequence (i.e. with degenerate bases)
#'   of the reverse primer.}
#'   \item{Degeneracy_rev}{Number of variants of the reverse primer.}
#' }
#'
#' If the option \code{showAllVariants == TRUE} was used for primer design
#' in \code{get_oligos()}, the following columns are also added:
#'
#' \describe{
#'   \item{All_fwd}{Lists with all sequence variants of the forward primer.}
#'   \item{GC_all_fwd}{Lists with the GC content of all
#'   sequence variants of the forward primer.}
#'   \item{Tm_all_fwd}{Lists with the Tm of all sequence variants of
#'   the reverse primer.}
#'   \item{All_rev}{Lists with all sequence variants of the reverse primer.}
#'   \item{GC_all_rev}{Lists with the GC content of all
#'   sequence variants of the forward primer.}
#'   \item{Tm_all_rev}{Lists with the Tm of all sequence variants of
#'   the reverse primer.}
#' }
#'
#' If a probe is used, the following columns are also included:
#'
#' \describe{
#'   \item{Begin_pr}{Position where the probe begins.}
#'   \item{End_pr}{Position where the probe ends.}
#'   \item{Length_pr}{Length of the probe.}
#'   \item{Majority_pr}{Majority sequence of the probe.}
#'   \item{GC_majority_pr}{GC-content of the probe
#'   (majority sequence), proportion.}
#'   \item{Tm_majority_pr}{Tm of the probe
#'   (majority sequence).}
#'   \item{Identity_pr}{Average identity of the probe.}
#'   \item{IUPAC_pr}{IUPAC sequence (i.e. with degenerate bases)
#'   of the probe.}
#'   \item{Degeneracy_pr}{Number of variants of the probe.}
#'   \item{Sense_pr}{Sense of the probe (pos or neg). If both probes are valid,
#'   the probe with the least G:s is selected.}
#' }
#'
#' If the option \code{showAllVariants == TRUE} was used for probe design
#' in \code{get_oligos()}, the following columns are also added:
#'
#' \describe{
#'   \item{All_pr}{Lists with all sequence variants of the probe.}
#'   \item{GC_all_pr}{Lists with the GC content of all
#'   sequence variants of the probe.}
#'   \item{Tm_all_pr}{Lists with the Tm of all sequence variants of
#'   the probe.}
#' }
#'
#' \describe{
#'   \item{Tm_difference_primer_probe}{Difference in Tm between the average
#'   Tm of the primer pair and the probe, majority sequences.}
#' }
#'
#' @seealso getOligos
#'
#' @export
getAssays <- function(primers,
                      probes = NULL,
                      length = 65:120,
                      maxTmDifferencePrimers = 2,
                      tmDifferenceProbes = NULL
                      ) {
  assays <- .combinePrimers(
    primers = primers, length = length,
    maxTmDifferencePrimers = maxTmDifferencePrimers)
  if (!is.null(probes)) {
    assays <- .addProbes(
      assays = assays, probes = probes, tmDifferenceProbes = tmDifferenceProbes
    )
  }
  assays
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
                            maxTmDifferencePrimers = 2
                            ) {
  if (!(maxTmDifferencePrimers > 0 && maxTmDifferencePrimers < 30)) {
    stop("'maxTmDifferencePrimers' must be from 0 to 30.", call. = FALSE)
  }
  if (!(min(length) >= 40 && max(length) <= 5000)) {
    stop("'length' must be from 40 to 5000.", call. = FALSE)
  }
  # Get all potential candidates for fwd and rev primers
  fwd <- primers[!is.na(primers$Majority), ]
  rev <- primers[!is.na(primers$Majority_RC), ]
  # Get all possible combinations of fwd and rev primers
  combinations <- expand.grid(
    fwd$Majority, rev$Majority_RC, stringsAsFactors = FALSE
  )
  names(combinations) <- c("fwdMajority", "revMajority")
  # Get indexes of the fwd and rev primer sequences
  indexFwd <- match(combinations$fwdMajority, primers$Majority)
  indexRev <- match(combinations$revMajority, primers$Majority_RC)
  # Make two datasets of x, one for fwd and one for rev
  fwd <- primers[indexFwd, ]
  rev <- primers[indexRev, ]
  # Add a tag on colnames before combining the two datasets
  colnames(fwd) <- paste0(colnames(fwd), "_fwd")
  colnames(rev) <- paste0(colnames(rev), "_rev")
  # Combine the two datasets
  assays <- dplyr::bind_cols(fwd, rev)
  assays <- tibble::as_tibble(assays)
  # Add amplicon length, tm difference and total degeneracy to the data
  Amplicon_length <- assays$End_rev - assays$Begin_fwd + 1
  Amplicon_length <- as.integer(Amplicon_length)
  Tm_difference_primer <- assays$Tm_majority_fwd - assays$Tm_majority_rev
  Tm_difference_primer <- abs(Tm_difference_primer)
  Begin <- assays$Begin_fwd
  End <- assays$End_rev
  Total_degeneracy <- assays$Degeneracy_fwd + assays$Degeneracy_rev
  Mean_identity <- mean(c(assays$Identity_fwd, assays$Identity_rev))
  assays <- tibble::add_column(
    assays, Begin, End, Amplicon_length,
    Tm_difference_primer, Mean_identity, Total_degeneracy, .before = "Begin_fwd"
  )
  # Drop columns that we do no longer need
  drop <- c("Majority_RC_fwd", "IUPAC_RC_fwd", "Majority_rev", "IUPAC_rev")
  assays <- assays[!(names(assays) %in% drop)]
  # Rename columns
  names(assays)[grep("_rc", names(assays))] <- c("Majority_rev", "IUPAC_rev")
  # Collect assays with desired amplicon length and tm difference
  assays <- assays[assays$Amplicon_length >= min(length), ]
  assays <- assays[assays$Amplicon_length <= max(length), ]
  assays <- assays[assays$Tm_difference_primer <= maxTmDifferencePrimers, ]
  if (nrow(assays) == 0L)
    stop("No assays were found.", call. = FALSE)
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
.addProbes <- function(assays, probes, tmDifferenceProbes = c(0, 20)) {
  if (is.null(tmDifferenceProbes)) tmDifferenceProbes <- c(-20, 20)
  if (!(min(tmDifferenceProbes) >= -20 && max(tmDifferenceProbes) <= 20)) {
    stop(
      "'tm_difference' must be from -20 to 20, e.g. c(-1, 5).",
      call. = FALSE
    )
  }
  # For each assay, take all probes that bind within the assay region
  probeCandidates <- purrr::map(seq_len(nrow(assays)), function(i) {
    # The probe has to begin after the fwd primer ends
    # and we want at least one base between (hence the + 2)
    from <- assays$End_fwd[[i]] + 2
    # The probe has to end before the rev primer begins
    # and we want at least one base in-between
    to <- assays$Begin_rev[[i]] - 2
    # Take all probes that begin and end 'within' the assay region
    probe <- probes[probes$Begin >= from & probes$End <= to, ]
    probe
  })
  # For each assay, we check how many probe candidates we have
  numberOfProbes <- purrr::map_int(probeCandidates, nrow)
  # Then, we need to pick all the assays that can harbor a probe.
  # For assays that has n number of probes, we need to
  # repeat that row n times.
  rowsToSelect <- purrr::map(
    seq_along(numberOfProbes), ~rep(.x, numberOfProbes[[.x]])
  )
  rowsToSelect <- unlist(rowsToSelect, use.names = FALSE)
  assays <- assays[rowsToSelect, ]
  if (nrow(assays) == 0L) {
    stop("No assays with probes could be generated.", call. = FALSE)
  }
  # Now, we can make a data frame of the probe candidates
  probeCandidates <- do.call("rbind", probeCandidates)
  # Select sense of the probe
  Sense <- purrr::map2_chr(
    probeCandidates$Majority, probeCandidates$Majority_RC, function(x, y) {
      if (is.na(x)) Sense <- "Neg"
      if (is.na(y)) Sense <- "Pos"
      if (!is.na(x) && !is.na(y)) {
        # If both probes are valid, select the probe with the least Gs
        gContentPos <- .gContent(x)
        gContentNeg <- .gContent(y)
        Sense <- ifelse(gContentPos <= gContentNeg, "Pos", "Neg")
      }
      Sense
  })
  probeCandidates$Majority <- ifelse(
    Sense == "Pos", probeCandidates$Majority, probeCandidates$Majority_RC
  )
  probeCandidates$IUPAC <- ifelse(
    Sense == "Pos", probeCandidates$IUPAC, probeCandidates$IUPAC_RC
  )
  probeCandidates <- tibble::add_column(probeCandidates, Sense)
  drop <- c("majority_rc", "iupac_rc")
  probeCandidates <- probeCandidates[!names(probeCandidates) %in% drop]
  names(probeCandidates) <- paste0(names(probeCandidates), "_pr")
  assays <- dplyr::bind_cols(assays, probeCandidates)
  assays$Mean_identity <- mean(
    c(assays$Identity_fwd, assays$Identity_rev, assays$Identity_pr)
  )
  assays$Total_degeneracy <- assays$Total_degeneracy + probeCandidates$Degeneracy_pr
  Tm_difference_primer_probe <- purrr::map_dbl(
    seq_len(nrow(assays)), function(x) {
      assays$Tm_majority_pr[[x]] - mean(
        assays$Tm_majority_fwd[[x]], assays$Tm_majority_rev[[x]]
      )
    })
  assays <- tibble::add_column(
    assays, Tm_difference_primer_probe, .after = "Tm_difference_primer"
  )
  assays <- assays[assays$Tm_difference_primer_probe >= min(tmDifferenceProbes), ]
  assays <- assays[assays$Tm_difference_primer_probe <= max(tmDifferenceProbes), ]
  if (nrow(assays) == 0L) {
    stop("No assays with probes could be generated.", call. = FALSE)
  }
  assays
}
