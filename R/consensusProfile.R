#' Get sequence information from an alignment
#'
#' \code{consensusProfile()} takes a DNA multiple alignment as input and
#' returns all the data needed for subsequent primer and probe design.
#' The function is a wrapper to
#' \code{Biostrings::consensusMatrix()} (Pages et al., 2020).
#'
#' @param x
#' A \code{Biostrings::DNAMultipleAlignment} object.
#'
#' @param ambiguityThreshold
#' "Detection level" for ambiguous bases.
#' All DNA bases that occur with a relative frequency higher than the
#' specified value will be included when the IUPAC consensus character
#' is determined.
#' Can range from `0 to 0.2, defaults to \code{0}.
#'
#' @return
#' An \code{RprimerProfile} object.
#'
#' @section Output:
#'
#' The output contains the following information:
#'
#' \describe{
#'   \item{position}{Position in the alignment.
#'   Note that masked columns in the original alignment are removed, and
#'   hence not taken into account when the position is determined.}
#'   \item{a}{Proportion of A.}
#'   \item{c}{Proportion of C.}
#'   \item{g}{Proportion of G.}
#'   \item{t}{Proportion of T.}
#'   \item{other}{Proportion of bases other than A, C, G, T.}
#'   \item{gaps}{Proportion of gaps (recognized as "-" in the alignment).}
#'   \item{majority}{Majority consensus sequence.
#'   Denotes the most frequently occurring nucleotide.
#'   If two or more bases occur with the same frequency,
#'   the consensus nucleotide will be randomly selected among these.}
#'   \item{identity}{Proportion of sequences, among all sequences with a
#'   DNA base (i.e., A, C, G or T), that has the majority consensus base.}
#'   \item{iupac}{The consensus sequence expressed in IUPAC format.
#'   The IUPAC consensus sequence only
#'   takes 'A', 'C', 'G', 'T' and '-' as input. Degenerate bases
#'   will be skipped. If a position only contains
#'   degenerate bases, the IUPAC consensus will be \code{NA} at that
#'   position.}
#'   \item{entropy}{Shannon entropy.
#'   Shannon entropy is a measurement of
#'   variability (Shannon, 1951). It is calculated among occurring DNA bases
#'   (gaps and ambiguous bases are not included) at each
#'   position in the alignment. A value of \code{0} indicate complete
#'   conservation and a high value indicate high variability.}
#'   \item{coverage}{Proportion of sequences in the target alignment,
#'   among all sequences with a DNA base, that are covered the IUPAC consensus
#'   character.
#'   The value will be 1 if there are no "remaining" DNA bases (and/or if
#'   \code{ambiguityThreshold = 0}).}
#' }
#'
#' @references
#'
#' Pages, H., Aboyoun, P., Gentleman R., and DebRoy S. (2020). Biostrings:
#' Efficient manipulation of biological strings. R package version
#' 2.57.2.
#'
#' Shannon, SE. (1951). Prediction and Entropy of Printed English.
#' Bell System Technical Journal 30 (1): 50-64.
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
    x <- .dnaBasesOnly(x)
    basesToInclude <- apply(x, 2, function(y) {
        paste(rownames(x)[y > ambiguityThreshold], collapse = ",")
    })
    basesToInclude <- unname(basesToInclude)
    consensus <- vapply(
        basesToInclude, .asIUPAC, character(1L),
        USE.NAMES = FALSE
    )
    if (any(is.na(consensus))) {
        warning(
            "The consensus sequence contains 'NA' at at least one position.
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
    x <- .dnaBasesOnly(x)
    identity <- apply(x, 2, max)
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
    x <- .dnaBasesOnly(x)
    entropy <- apply(x, 2, function(y) ifelse(y == 0, 0, y * log2(y)))
    entropy <- abs(colSums(entropy))
    entropy[is.na(entropy)] <- 0
    unname(entropy)
}

#' Coverage
#'
#' \code{.coverage()} calculates the proportion of bases
#' that are covered within the ambiguous (IUPAC) bases. Gaps and ambiguous
#' bases are not included in the calculation.
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
    x <- .dnaBasesOnly(x)
    x[x > ambiguityThreshold] <- 0
    coverage <- 1 - colSums(x)
    coverage[is.na(coverage)] <- 1
    unname(coverage)
}
