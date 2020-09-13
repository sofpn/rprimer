#' Count the degeneracy of a DNA sequence
#'
#' \code{count_degenerates} counts the number of unique sequences of
#' a DNA sequence with degenerate bases.
#'
#' @param x a DNA sequence (a character vector of length one, e.g. 'cttgg').
#'
#' @details
#' Valid bases for \code{x} are 'a', 'c', 'g', 't', 'r', 'y', 'm',
#' 'k', 's', 'w', n', 'h', 'd', 'v', 'b' and '-'.
#'
#' @return The number of unique sequences of x (an integer).
#'
#' @examples
#' count_degeneracy("cttnra")
#'
#' @export
count_degeneracy <- function(x) {
  if (typeof(x) != "character") {
    stop("'x' must be a character vector.", call. = FALSE)
  }
  if (grepl("[^acgtrymkswnhdvb-]", x)) {
    stop("'x' contains at least one invalid base. \n
      Valid bases are 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w',
      'n', 'h', 'd', 'v', 'b' and '-'",
         call. = FALSE
    )
  }
  x <- split_sequence(x)
  # Find the number of nucleotides at each position in x
  n_nucleotides <- degeneracy_lookup[x]
  # Calculate the total number of DNA sequences in x
  degeneracy <- prod(n_nucleotides)
  degeneracy
}
