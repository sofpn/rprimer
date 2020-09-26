#' Count the degeneracy of a DNA sequence
#'
#' \code{countDegeneracy} returns the number of unique variants of
#' a DNA sequence with degenerate bases.
#'
#' @param x
#' A DNA sequence (a character vector of length one).
#' Valid bases are ACGTRYSWKMBDHVN-.
#'
#' @return The number of unique sequences of x (an integer).
#'
#' @examples
#' countDegeneracy("GTTTCCRT")
#'
#' @keywords internal
countDegeneracy <- function(x) {
  if (typeof(x) != "character" || length(x) != 1) {
    stop("'x' must be a character vector of length one.", call. = FALSE)
  }
  x <- toupper(x)
  if (grepl(paste0("[^", allBases, "]"), x)) {
    stop(paste0("'x' can only contain bases ", allBases, "."), call. = FALSE)
  }
  x <- splitSequence(x)
  nNucleotides <- degeneracyLookup[x]
  degeneracy <- prod(nNucleotides)
  degeneracy
}
