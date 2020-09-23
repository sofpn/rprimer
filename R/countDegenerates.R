#' Count the number of degenerate bases in a DNA sequence
#'
#' \code{count_degenerates} returns the number of degenerate bases in a DNA
#' sequence.
#'
#' @param x a DNA sequence (a character vector of length one, e.g. 'cttgg').
#'
#' @details Valid bases for \code{x} are 'a', 'c', 'g', 't', 'r', 'y', 'm',
#' 'k', 's', 'w', n', 'h', 'd', 'v', 'b' and '-'.
#'
#' @return the number of degenerate bases in \code{x} (an integer).
#'
#' @examples
#' count_degenerates("cttnra")
#'
#' @export
count_degenerates <- function(x) {
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
  nt <- c("a", "c", "g", "t", "-")
  x <- split_sequence(x)
  count <- length(x[!x %in% nt])
  count
}
