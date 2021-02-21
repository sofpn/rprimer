#' Get sequence information from an alignment
#'
#' \code{consensusProfile()} takes a DNA multiple alignment as input and
#' returns all the data needed for subsequent primer and probe design.
#'
#' @param x
#' A \code{Biostrings::DNAMultipleAlignment} object.
#'
#' @param ambiguityThreshold
#' A number [0, 0.2], defaults to 0.
#' At each position, all nucleotides with a proportion
#' higher than the \code{ambiguityThreshold} will be included in
#' the IUPAC consensus character.
#'
#' @return
#' An \code{RprimerProfile} object, which contains the following information:
#'
#' \describe{
#'   \item{position}{Position in the alignment. Note that masked columns
#'   from the original alignment are removed, and
#'   hence not taken into account when position is determined.}
#'   \item{a}{Proportion of A.}
#'   \item{c}{Proportion of C.}
#'   \item{g}{Proportion of G.}
#'   \item{t}{Proportion of T.}
#'   \item{other}{Proportion of bases other than A, C, G, T.}
#'   \item{gaps}{Proportion of gaps (recognized as "-" in the alignment).}
#'   \item{majority}{Majority consensus sequence.
#'   The most frequently occurring nucleotide.
#'   If two or more bases occur with the same frequency,
#'   the consensus nucleotide will be randomly selected among these bases.}
#'   \item{identity}{Proportion of the most common nucleotide.
#'   Gaps (-), as well as bases other than A, C, G and T are excluded from the
#'   calculation.}
#'   \item{iupac}{
#'   The consensus sequence expressed in IUPAC format.
#'   Note that the IUPAC consensus sequence only
#'   takes 'A', 'C', 'G', 'T' and '-' as input. Degenerate bases
#'   present in the alignment will be skipped. If a position only contains
#'   degenerate/invalid bases, the IUPAC consensus will be \code{NA} at that
#'   position.}
#'   \item{entropy}{Shannon entropy.
#'   Shannon entropy is a measurement of
#'   variability.
#'   First, for each nucleotide that occurs at a specific position,
#'   \code{p*log2(p)}, is calculated, where \code{p} is the proportion of
#'   that nucleotide. Then, these values are summarized,
#'   and multiplied by \code{-1}.
#'   A value of \code{0} indicate no variability and a high value
#'   indicate high variability.
#'   Gaps (-), as well as bases other than
#'   A, C, G and T are excluded from the calculation.}
#'   \item{coverage}{The proportion of bases that are included the
#'   consensus/ambiguous (IUPAC) base.
#'   Will be one if there are no "remaining" bases (and if
#'   \code{ambiguityThreshold = 0}).
#'   Gaps (-), as well as bases other than A, C, G and T are excluded from the
#'   calculation.}
#' }
#'
#' @references
#' This function is a wrapper to \code{Biostrings::consensusMatrix()}:
#'
#' H. Pages, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings:
#' Efficient manipulation of biological strings. R package version
#' 2.57.2.
#'
#' @export
#'
#' @examples
#' data("exampleRprimerAlignment")
#' consensusProfile(exampleRprimerAlignment)
consensusProfile <- function(x, ambiguityThreshold = 0) {
    if (!methods::is(x, "DNAMultipleAlignment")) {
        stop("'x' must be a DNAMultipleAlignment object.", call. = FALSE)
    }
    if (!(ambiguityThreshold >= 0 && ambiguityThreshold <= 0.2)) {
        stop(
            paste0("'ambiguityThreshold' must be a number from 0 to 0.2."),
            call. = FALSE
        )
    }
    x <- .consensusMatrix(x)
    profile <- list()
    profile$position <- seq_len(ncol(x))
    profile$a <- unname(x["A", ])
    profile$c <- unname(x["C", ])
    profile$g <- unname(x["G", ])
    profile$t <- unname(x["T", ])
    profile$other <- unname(x["other", ])
    profile$gaps <- unname(x["-", ])
    profile$majority <- .majorityConsensus(x)
    profile$identity <- .nucleotideIdentity(x)
    profile$iupac <- .iupacConsensus(x, ambiguityThreshold)
    profile$iupac[profile$majority == "-"] <- "-"
    profile$entropy <- .shannonEntropy(x)
    profile$coverage <- .coverage(x, ambiguityThreshold)
    RprimerProfile(profile)
}

# Helpers ======================================================================

#' Get consensus matrix
#'
#' @param x A \code{Biostrings::DNAMultipleAlignment} object.
#'
#' @return A consensus matrix.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerAlignment")
#' .consensusMatrix(exampleRprimerAlignment)
.consensusMatrix <- function(x) {
    x <- Biostrings::consensusMatrix(x, as.prob = TRUE)
    x <- x[, colSums(!is.na(x)) > 0, drop = FALSE] ## Removes masked columns
    colnames(x) <- seq_len(ncol(x))
    x <- x[(rownames(x) != "+" & rownames(x) != "."), , drop = FALSE]
    bases <- c("A", "C", "G", "T", "-")
    other <- colSums(x[!rownames(x) %in% bases, , drop = FALSE])
    x <- x[rownames(x) %in% bases, , drop = FALSE]
    rbind(x, other)
}

#' Find the most common base in a named vector
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerAlignment")
#' x <- .consensusMatrix(exampleRprimerAlignment)
#' .findMostCommonBase(x[, 1])
.findMostCommonBase <- function(x) {
    mostCommon <- names(x)[x == max(x)]
    if (length(mostCommon > 1)) mostCommon <- sample(mostCommon, 1)
    mostCommon
}

#' Majority consensus sequence
#'
#' @param x A consensus matrix.
#'
#' @return The majority consensus sequence (a character vector).
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerAlignment")
#' x <- .consensusMatrix(exampleRprimerAlignment)
#' .majorityConsensus(x)
.majorityConsensus <- function(x) {
    consensus <- apply(x, 2, .findMostCommonBase)
    unname(consensus)
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
#' \code{.asIUPAC()} will return NA.
#'
#' @return The corresponding IUPAC base.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' .asIUPAC("A,G,C")
.asIUPAC <- function(x) {
    x <- unlist(strsplit(x, split = ","), use.names = FALSE)
    x <- x[order(x)]
    x <- unique(x)
    bases <- c("A", "C", "G", "T", "-")
    x <- x[x %in% bases]
    x <- paste(x, collapse = ",")
    unname(lookup$iupac[x])
}

#' Subset a consensus matrix
#'
#' @param x A consensus matrix.
#'
#' @return A consensus bases with DNA bases (A, C, G, T) only.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerAlignment")
#' x <- .consensusMatrix(exampleRprimerAlignment)
#' .dnaBasesOnly(x)
.dnaBasesOnly <- function(x) {
    bases <- c("A", "C", "G", "T")
    s <- x[rownames(x) %in% bases, , drop = FALSE]
    apply(s, 2, function(x) x / sum(x))
}

#' IUPAC consensus sequence
#'
#' @param x A consensus matrix.
#'
#' @param ambiguityThreshold
#' At each position, all nucleotides with a proportion
#' higher than the threshold will be included in
#' the IUPAC consensus sequence.
#'
#' @return The consensus sequence (a character vector).
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerAlignment")
#' x <- .consensusMatrix(exampleRprimerAlignment)
#' .iupacConsensus(x)
.iupacConsensus <- function(x, ambiguityThreshold = 0) {
    s <- .dnaBasesOnly(x)
    basesToInclude <- apply(s, 2, function(x) {
        paste(rownames(s)[x > ambiguityThreshold], collapse = ",")
    })
    basesToInclude <- unname(basesToInclude)
    consensus <- vapply(
        basesToInclude, .asIUPAC, character(1),
        USE.NAMES = FALSE
    )
    if (any(is.na(consensus))) {
        warning(
            "The consensus sequence contains 'NA' at at least one positions.
        Perharps an ambiguos base was the most common base
        at these positions.",
            call. = FALSE
        )
    }
    consensus
}

#' Nucleotide identity
#'
#' @param x A consensus matrix.
#'
#' @return The nucleotide identity (a numeric vector).
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerAlignment")
#' x <- .consensusMatrix(exampleRprimerAlignment)
#' .nucleotideIdentity(x)
.nucleotideIdentity <- function(x) {
    s <- .dnaBasesOnly(x)
    identity <- apply(s, 2, max)
    identity <- unname(identity)
    identity[is.na(identity)] <- 1
    identity
}

#' Shannon entropy
#'
#' @param x A consensus matrix.
#'
#' @return The Shannon entropy (a numeric vector).
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerAlignment")
#' x <- .consensusMatrix(exampleRprimerAlignment)
#' .shannonEntropy(x)
.shannonEntropy <- function(x) {
    s <- .dnaBasesOnly(x)
    entropy <- apply(s, 2, function(x) ifelse(x == 0, 0, x * log2(x)))
    entropy <- abs(colSums(entropy))
    entropy <- unname(entropy)
    entropy[is.na(entropy)] <- 0
    entropy
}

#' Coverage
#'
#' \code{.coverage()} calculates the proportion of bases
#' that are covered within the ambiguous (IUPAC) bases. Gaps as well as bases
#' other than A, C, G and T are not included in the calculation.
#'
#' @param x A consensus matrix.
#'
#' @return The coverage (a numeric vector). A value of 1 means that
#' all bases are covered within the ambiguous base.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerAlignment")
#' x <- .consensusMatrix(exampleRprimerAlignment)
#' .coverage(x, ambiguityThreshold = 0.05)
.coverage <- function(x, ambiguityThreshold = 0) {
    s <- .dnaBasesOnly(x)
    s[s > ambiguityThreshold] <- 0
    coverage <- 1 - colSums(s)
    coverage[is.na(coverage)] <- 1
    unname(coverage)
}
