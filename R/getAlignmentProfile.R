## How to deal with roi from iranges?

#' Get the sequence profile of an alignment
#'
#' \code{getAlignmentProfile()} returns a matrix with the
#' proportion of each nucleotide at each position within an alignment
#' of DNA sequences. The function is a wrapper around
#' \code{Biostrings::consensusMatrix()}.
#'
#' @param x
#' A \code{Biostrings::DNAMultipleAlignment} object.
#'
#' @return
#' An \code{RprimerProfile} object, containing a numeric matrix with the
#' proportion of each
#' nucleotide at each position within an alignment
#' of DNA sequences. The matrix has six rows,
#' named 'A', 'C', 'G', 'T', '-' and 'Other'. '-' represents gaps and
#' 'Other' represents nucleotides other than A, C, G and T, e.g.
#' wobble bases. The columns are named
#' according to which position they correspond to in the alignment.
#' To do: Note that
#' columns with NA are removed.
#'
#'
#' @examples
#' data("exampleRprimerAlignment")
#' getAlignmentProfile(exampleRprimerAlignment)
#'
#' @references
#' H. Pagès, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings:
#' Efficient manipulation of biological strings. R package version
#' 2.57.2.
#'
#' @export
getAlignmentProfile <- function(x) {
    if (!methods::is(x, "DNAMultipleAlignment")) {
        stop("'x' must be a DNAMultipleAlignment object.", call. = FALSE)
    }
    x <- Biostrings::consensusMatrix(x, as.prob = TRUE)
    x <- x[, colSums(!is.na(x)) > 0]
    colnames(x) <- seq_len(ncol(x))
    x <- x[(rownames(x) != "+" & rownames(x) != "."), ]
    bases <- c("A", "C", "G", "T", "-")
    other <- colSums(x[!rownames(x) %in% bases, ])
    x <- x[rownames(x) %in% bases, ]
    x <- rbind(x, other)
    x <- RprimerProfile(x)
    x
}
