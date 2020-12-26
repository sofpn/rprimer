#' Get sequence information from an alignment
#'
#' \code{getConsensusProfile()} returns sequence information from an alignment
#' of DNA sequences.
#'
#' @param x
#' A \code{Biostrings::DNAMultipleAlignment} object.
#'
#' @param iupacThreshold
#' A number [0, 0.2].
#' At each position, all nucleotides with a proportion
#' higher than the \code{iupacThreshold} will be included in
#' the IUPAC consensus sequence. Defaults to 0.
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
#'   \item{gaps}{proportion of gaps (recognized as "-" in the alignment).}
#'   \item{majority}{Majority consensus sequence.
#'   The most frequently occurring nucleotide.
#'   If two or more bases occur with the same frequency,
#'   the consensus nucleotide will be randomly selected among these bases.}
#'   \item{identity}{Proportion of the most common nucleotide.
#'   Gaps (-), as well as bases other than A, C, G and T are excluded from the
#'   calculation.}
#'   \item{iupac}{
#'   The consensus sequence expressed in IUPAC format (i.e. with wobble bases)
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
#' }
#'
#' @examples
#' data("exampleRprimerAlignment")
#' getConsensusProfile(exampleRprimerAlignment)
#'
#' @references
#' This function is a wrapper around \code{Biostrings::consensusMatrix()}:
#'
#' H. Pagès, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings:
#' Efficient manipulation of biological strings. R package version
#' 2.57.2.
#'
#' @export
getConsensusProfile <- function(x, iupacThreshold = 0) {
    if (!methods::is(x, "DNAMultipleAlignment")) {
        stop("'x' must be a DNAMultipleAlignment object.")
    }
    x <- .getConsensusMatrix(x)
    position <- seq_len(ncol(x))
    a <- unname(x["A", ])
    c <- unname(x["C", ])
    g <- unname(x["G", ])
    t <- unname(x["T", ])
    other <- unname(x["other", ])
    gaps <- unname(x["-", ])
    majority <- .majorityConsensus(x)
    identity <- .nucleotideIdentity(x)
    iupac <- .iupacConsensus(x, iupacThreshold = iupacThreshold)
    iupac[majority == "-"] <- "-"
    entropy <- .shannonEntropy(x)
    df <- data.frame(
        position, a, c, g, t, other, gaps, majority, identity, iupac, entropy
    )
    RprimerProfile(df)
}

# Helpers =====================================================================

#' Get consensus matrix
#'
#' @param x A \code{Biostrings::DNAMultipleAlignment} object.
#'
#' @return A consensus matrix.
#'
#' @keywords internal
#'
#' @noRd
.getConsensusMatrix <- function(x) {
    x <- Biostrings::consensusMatrix(x, as.prob = TRUE)
    x <- x[, colSums(!is.na(x)) > 0] # Removes the masked columns
    colnames(x) <- seq_len(ncol(x))
    x <- x[(rownames(x) != "+" & rownames(x) != "."), ]
    bases <- c("A", "C", "G", "T", "-")
    other <- colSums(x[!rownames(x) %in% bases, ])
    x <- x[rownames(x) %in% bases, ]
    rbind(x, other)
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
.majorityConsensus <- function(x) {
    .findMostCommonBase <- function(x, y) {
        mostCommon <- rownames(x)[y == max(y)]
        if (length(mostCommon > 1)) mostCommon <- sample(mostCommon, 1)
        mostCommon
    }
    consensus <- apply(x, 2, function(y) .findMostCommonBase(x, y))
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
    bases <- c("A", "C", "G", "T", "-")
    x <- x[x %in% bases]
    x <- paste(x, collapse = ",")
    unname(lookup$iupac[x])
}

#' IUPAC consensus sequence
#'
#' @param x A consensus matrix.
#'
#' @param iupacThreshold
#' Optional. A number [0, 0.2]
#' At each position, all nucleotides with a proportion
#' higher than the threshold will be included in
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
.iupacConsensus <- function(x, iupacThreshold = 0) {
    if (
        !is.numeric(iupacThreshold) ||
            iupacThreshold < 0 || iupacThreshold > 0.2
    ) {
        stop(
            paste0("'iupacThreshold' must be a number from 0 to 0.2."),
            call. = FALSE
        )
    }
    bases <- c("A", "C", "G", "T", "-")
    x <- x[rownames(x) %in% bases, ]
    basesToInclude <- apply(x, 2, function(y) {
        paste(rownames(x)[y > iupacThreshold], collapse = ",")
    })
    basesToInclude <- unname(basesToInclude)
    consensus <- purrr::map_chr(basesToInclude, ~ .asIUPAC(.x))
    if (any(is.na(consensus))) {
        warning("The consensus sequence contain NAs. \n
    Try to lower the 'iupacThreshold' value.", call. = FALSE)
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
.nucleotideIdentity <- function(x) {
    bases <- c("A", "C", "G", "T")
    s <- x[rownames(x) %in% bases, ]
    s <- apply(s, 2, function(x) x / sum(x))
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
.shannonEntropy <- function(x) {
    bases <- c("A", "C", "G", "T")
    s <- x[rownames(x) %in% bases, ]
    s <- apply(s, 2, function(x) x / sum(x))
    entropy <- apply(s, 2, function(x) ifelse(x == 0, 0, x * log2(x)))
    entropy <- -colSums(entropy)
    entropy <- unname(entropy)
    entropy <- abs(entropy)
    entropy[is.na(entropy)] <- 0
    entropy
}
