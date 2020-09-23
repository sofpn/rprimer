#' Count the number of degenerate bases in a DNA sequence
#'
#' \code{countDegenerates} returns the number of degenerate bases in a DNA
#' sequence.
#'
#' @param x
#' A DNA sequence (a character vector of length one).
#' Valid bases are ACGTRYSWKMBDHVN-.
#'
#' @return The number of degenerate bases in \code{x} (an integer).
#'
#' @examples
#' count_degenerates("CTTRNA")
#'
#' @export
countDegenerates <- function(x) {
  if (typeof(x) != "character" || length(x) != 1) {
    stop("'x' must be a character vector of length one.", call. = FALSE)
  }
  x <- toupper(x)
  if (grepl(paste0("[^", allBases, "]"), x)) {
    stop(paste0("'x' can only contain bases ", allBases, "."), call. = FALSE)
  }
  nt <- c("A", "C", "G", "T", "-")
  x <- splitSequence(x)
  count <- length(x[!x %in% nt])
  count
}
