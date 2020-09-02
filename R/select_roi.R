#' Select region of interest
#'
#' \code{select_roi} selects a specified region of interest within
#' an alignment of DNA sequences.
#'
#' @param  x
#' An alignment of DNA sequences (an object of class 'rprimer_alignment').
#'
#' @param from Where the roi begins (an integer). The default is 1.
#'
#' @param to
#' Where the roi ends (an integer). The default is \code{NULL}.
#' In that case, \code{to} will be set as the last position in the alignment.
#'
#' @return
#' The roi of the alignment
#' (an object of class 'rprimer_alignment').
#'
#' @examples
#' # Select the first 1000 bases
#' select_roi(example_rprimer_alignment, from = 1, to = 1000)
#'
#' @export
select_roi <- function(x, from = 1, to = NULL) {
  if (!is.rprimer_alignment(x)) {
    stop("'x' must be an rprimer_alignment object.", call. = FALSE)
  }
  splitted <- purrr::map(x, split_sequence)
  # Get the length of the alignment
  # All sequences are of equal length in an rprimer_alignment object
  aln_length <- length(splitted[[1]])
  if (is.null(to)) to <- aln_length
  if (!(is.numeric(from) && is.numeric(to))) {
    stop(paste0(
      "'from' and 'to' must be positive integers. \n
      You've set 'from' to ", from, " and 'to' to ", to, "."), call. = FALSE
    )
  }
  if (to > aln_length) {
    stop(
      "You've set 'to' at a position longer than the alignment length.",
      call. = FALSE
    )
  }
  if (to > aln_length) {
    stop(
      "You've set 'from' at a position longer than the alignment length.",
      call. = FALSE
    )
  }
  if (to < from) {
    stop("'to' must be greater than 'from'", call. = FALSE)
  }
  if (!(from >= 1 && to >= 1)) {
    stop(paste0(
      "'from' and 'to' must be positive integers. \n
      You've set 'from' to ", from, " and 'to' to ", to, "."), call. = FALSE
    )
  }
  splitted <- purrr::map(splitted, ~.x[from:to])
  x <- purrr::map(splitted, paste, collapse = "")
  x <- new_rprimer_alignment(x)
  x
}
