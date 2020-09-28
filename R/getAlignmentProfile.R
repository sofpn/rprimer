#' Get the sequence profile of an alignment
#'
#' \code{getAlignmentProfile} returns a matrix with the
#' proportion of each nucleotide at each position within an alignment
#' of DNA sequences.The function is a wrapper around
#' Biostrings::consensusMatrix.
#'
#' @param x
#' An alignment of DNA sequences.
#'
#' @return
#' The sequence profile (an object of class RprimerProfile).
#' A numeric m x n matrix,
#' where m is the number of unique bases in the alignment, and n is the
#' number of positions in the alignment.
#'
#' @references
#' H. Pagès, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings:
#' Efficient manipulation of biological strings. R package version 2.57.2.
#'
#' @examples
#' data("exampleRprimerAlignment")
#' getAlignmentProfile(exampleRprimerAlignment)
#' @export
getAlignmentProfile <- function(x) {
    if (!methods::is(x, "DNAMultipleAlignment")) {
        stop("'x' must be a DNAMultipleAlignment object.", call. = FALSE)
    }
    x <- Biostrings::consensusMatrix(x, as.prob = TRUE)
    x <- x[, colSums(!is.na(x)) > 0]
    x <- x[(rownames(x) != "+" & rownames(x) != "."), ]
    bases <- c("A", "C", "G", "T", "-")
    other <- colSums(x[!rownames(x) %in% bases, ])
    x <- x[rownames(x) %in% bases, ]
    x <- rbind(x, other)
    colnames(x) <- seq_len(ncol(x))
    x <- RprimerProfile(x)
    x
}
