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
#' @examples
#'
#'
#' @export
getAlignmentProfile <- function(x) {
  # if (class(myAlignment) != "DNAMultipleAlignment")
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
