## To do: make S4 class

#' Get alignment properties
#'
#' \code{getAlignmentProperties()} returns sequence information from
#' an RprimerProfile object
#'
#' @param x An \code{RprimerProfile} object.
#'
#' @param iupacThreshold
#' A number (0, 0.2].
#' At each position, all nucleotides with a proportion
#' higher or equal to the \code{iupacThreshold} will be included in
#' the IUPAC consensus sequence. It defaults to 0.
#'
#' @section
#' Majority consensus sequence:
#' The most frequently occurring nucleotide.
#' If two or more bases occur with the same frequency,
#' the consensus nucleotide will be randomly selected among these bases.
#'
#' @section
#' IUPAC consensus sequence:
#' The consensus sequence expressed in IUPAC format (i.e. with wobble bases)
#' Note that the IUPAC consensus sequence only
#' takes 'A', 'C', 'G', 'T' and '-' as input. Degenerate bases
#' present in the alignment will be skipped. If a position only contains
#' degenerate/invalid bases, the IUPAC consensus will be \code{NA} at that
#' position.
#'
#' @section Gaps:
#' Proportion of gaps. Gaps are recognized as "-".
#'
#' @section Identity:
#' Proportion of
#' the most common base. Gaps (-),
#' as well as bases other than A, C, G and T are excluded from the
#' calculation.
#'
#' @section Entropy:
#' Shannon entropy is a measurement of
#' variability. First, for each nucleotide that occurs at a specific position,
#' \code{p*log2(p)}, is calculated, where \code{p} is the proportion of
#' that nucleotide. Then, these values are summarized,
#' and multiplied by \code{-1}.
#' A value of \code{0} indicate no variability and a high value
#' indicate high variability.
#' Gaps (-), as well as bases other than
#' A, C, G and T are excluded from the calculation.
#'
#' @return
#' A tibble (a data frame) with the following information:
#'
#' #' \describe{
#'   \item{Position}{Position in the alignment.}
#'   \item{Majority}{Majority consensus sequence.}
#'   \item{IUPAC}{IUPAC consensus sequence.}
#'   \item{Gaps}{Proportion of gaps.}
#'   \item{Identity}{Proportion of the most common nucleotide.}
#'   \item{Entropy}{Shannon entropy}
#' }
#'
#' @examples
#' data("exampleRprimerProfile")
#' getAlignmentProperties(exampleRprimerProfile)
#' getAlignmentProperties(exampleRprimerProfile, iupacThreshold = 0.1)
#'
#' @export
getAlignmentProperties <- function(x, iupacThreshold = 0) {
    if (!methods::is(x, "RprimerProfile")) {
        stop("'x' must be an RprimerProfile object", call. = FALSE)
    }
    x <- SummarizedExperiment::assay(x)
    position <- seq_len(ncol(x))
    majority <- .majorityConsensus(x)
    IUPAC <- .iupacConsensus(x, threshold = iupacThreshold)
    gaps <- .gapFrequency(x)
    identity <- .nucleotideIdentity(x)
    entropy <- .shannonEntropy(x)
    properties <- tibble::tibble(
        "Position" = position,
        "Majority" = majority,
        "IUPAC" = IUPAC,
        "Gaps" = gaps,
        "Identity" = identity,
        "Entropy" = entropy
    )
    properties
}

# Helpers =====================================================================

#' Majority consensus sequence
#'
#' \code{.majorityConsensus()} returns the majority consensus sequence of an
#' alignment of DNA sequences.
#'
#' @param x A numeric matrix.
#'
#' @return The majority consensus sequence (a character vector).
#'
#' @keywords internal
#'
#' @noRd
.majorityConsensus <- function(x) {
    .findMostCommonBase <- function(x, y) {
        mostCommon <- rownames(x)[y == max(y)]
        if (length(mostCommon > 1)) {
            mostCommon <- sample(mostCommon, 1)
        }
        mostCommon
    }
    consensus <- apply(x, 2, function(y) .findMostCommonBase(x, y))
    consensus <- unname(consensus)
    consensus
}

#' Convert DNA nucleotides into the corresponding IUPAC base
#'
#' \code{.asIUPAC()} takes several DNA nucleotides as input,
#' and returns the degenerate base in IUPAC format.
#'
#' @param x
#' A character vector of length one, containing DNA
#' nucleotides (valid bases are A, C, G, T, - and .). Each base must
#' be separated by a comma (,), e.g. 'A,C,G'.
#' Characters other than A, C, G, T, and - will be ignored. However, when
#' the input only consist of invalid bases,
#' or if the bases are not separated by ',',
#' \code{.asIUPAC} will return NA.
#'
#' @return The corresponding IUPAC base.
#'
#' @keywords internal
#'
#' @noRd
.asIUPAC <- function(x) {
    x <- toupper(x)
    x <- gsub(" ", "", x)
    x <- unlist(strsplit(x, split = ","), use.names = FALSE)
    x <- x[order(x)]
    x <- unique(x)
    bases <- c("A", "C", "G", "T", "-", ".")
    match <- x %in% bases
    x <- x[!(match == FALSE)]
    x <- paste(x, collapse = ",")
    iupacBase <- unname(rprimerGlobals$iupacLookup[x])
    iupacBase
}

#' IUPAC consensus sequence
#'
#' \code{.iupacConsensus()} returns the IUPAC consensus sequence from an
#' PrprimerProfile object.
#'
#' @param x A numeric matrix.
#'
#' @param threshold
#' Optional. A number (0, 0.2]
#' At each position, all nucleotides with a proportion
#' higher or equal to the threshold will be included in
#' the IUPAC consensus sequence. The default is 0.
#'
#' @return The consensus sequence (a character vector).
#'
#' @seealso \code{.asIUPAC} for further info on how the IUPAC consensus
#' sequence is determined.
#'
#' @keywords internal
#'
#' @noRd
.iupacConsensus <- function(x, threshold = 0) {
    if (!is.double(threshold) || threshold < 0 || threshold > 0.2) {
        stop(paste0(
            "'treshold' must be higher than 0 and less or equal to 0.2."
        ), call. = FALSE)
    }
    bases <- c("A", "C", "G", "T", "-")
    x <- x[rownames(x) %in% bases, ]
    basesToInclude <- apply(x, 2, function(y) {
        paste(rownames(x)[y > threshold], collapse = ",")
    })
    basesToInclude <- unname(basesToInclude)
    consensus <- purrr::map_chr(basesToInclude, ~.asIUPAC(.x))
    if (any(is.na(consensus))) {
        warning("The consensus sequence contain NAs. \n
    Try to lower the threshold value.", call. = FALSE)
    }
    consensus
}

#' Gap frequency
#'
#' @param x A numeric matrix.
#'
#' @return The gap frequency (a numeric vector).
#'
#' @keywords internal
#'
#' @noRd
.gapFrequency <- function(x) {
    if ("-" %in% rownames(x)) {
        gaps <- x[rownames(x) == "-", ]
        gaps <- unname(gaps)
    } else {
        gaps <- rep(0, ncol(x))
    }
    gaps
}

#' Nucleotide identity
#'
#' @param x A numeric matrix.
#'
#' @return The nucleotide identity (a numeric vector).
#' The nucleotide identity can range between (0, 1].
#'
#' @keywords internal
#'
#' @noRd
.nucleotideIdentity <- function(x) {
    bases <- c("A", "C", "G", "T")
    s <- x[rownames(x) %in% bases, ]
    s <- apply(s, 2, function(x) x / sum(x))
    identity <- apply(s, 2, max)
    identity <- unname(identity)
    identity[is.na(identity)] <- 0
    identity
}

#' Shannon entropy
#'
#' @param x A numeric matrix.
#'
#' @return The Shannon entropy (a numeric vector).
#'
#' @keywords internal
#'
#' @noRd
.shannonEntropy <- function(x) {
    bases <- c("A", "C", "G", "T")
    s <- x[rownames(x) %in% bases, ]
    s <- apply(s, 2, function(x) x / sum(x))
    entropy <- apply(s, 2, function(x) {
        ifelse(x == 0, 0, x * log2(x))
    })
    entropy <- -colSums(entropy)
    entropy <- unname(entropy)
    entropy <- abs(entropy)
    entropy[is.na(entropy)] <- 0
    entropy
}
