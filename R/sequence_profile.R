#' Get the sequence profile of an alignment
#'
#' \code{sequence_profile} returns a matrix with the
#' proportion of each nucleotide at each position within an alignment
#' of DNA sequences.
#'
#' @param x
#' An alignment of DNA sequences (an object of class 'rprimer_alignment').
#'
#' @return
#' The sequence profile (an object of class rprimer_profile').
#' A numeric m x n matrix,
#' where m is the number of unique bases in the alignment, and n is the
#' number of positions in the alignment.
#'
#' @examples
#' sequence_profile(example_rprimer_alignment)
#'
#' @export
sequence_profile <- function(x) {
  if (!is.rprimer_alignment(x)) {
    stop("'x' must be an rprimer_alignment object.", call. = FALSE)
  }
  splitted <- purrr::map(x, split_sequence)
  # Make a matrix
  matr <- do.call("rbind", splitted)
  # Get all unique bases in the dataset and sort them in alphabetical order
  bases <- unique(sort(unlist(apply(matr, 1, unique), use.names = FALSE)))
  # Count the occurence of each base at each position
  count_base <- function(x, base) length(x[which(x == base)])
  counts <- apply(matr, 2, function(x) {
    purrr::map_int(bases, ~count_base(x, base = .x))
  })
  # Present the data as proportions instead of counts
  proportions <- counts / nrow(matr)
  # Base as rowname
  rownames(proportions) <- bases
  # Nucleotide positon as colname
  colnames(proportions) <- seq_len(dim(matr)[[2]])
  proportions <- new_rprimer_profile(proportions)
  proportions
}
