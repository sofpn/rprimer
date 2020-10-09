#' Get the sequence profile of an alignment
#'
#' \code{getAlignmentProfile} returns a matrix with the
#' proportion of each nucleotide at each position within an alignment
#' of DNA sequences.The function is a wrapper around
#' Biostrings::consensusMatrix.
#'
#' @param x
#' A Biostrings::DNAMultipleAlignment object.
#'
#' @return
#' An RprimerProfile object, containing a numeric matrix with the
#' proportion of each
#' nucleotide at each position within an alignment
#' of DNA sequences. The matrix has six rows,
#' named 'A', 'C', 'G', 'T', '-' and 'Other'. '-' represents gaps and
#' 'Other' represents nucleotides other than A, C, G and T, e.g.
#' wobble bases. The columns are named
#' according to which position they correspond to in the input alignment.
#'
#'
#' @examples
#' data("exampleRprimerAlignment")
#' getAlignmentProfile(exampleRprimerAlignment)
#'
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
