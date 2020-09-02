#' Remove positions with high gap frequency
#'
#' \code{remove_gaps} removes all positions with a gap frequency
#' higher than a stated threshold in an alignment of DNA sequences.
#'
#' @param x
#' An alignment of DNA sequences (an object of class 'rprimer_alignment').
#'
#' @param threshold A number between 0.5 and <1 (the default is 0.5).
#'
#' @return
#' An alignment (an object class 'rprimer_alignment'), where
#' all positions with a gap proportion higher than the stated threshold
#' are removed.
#'
#' @examples
#' remove_gaps(example_rprimer_alignment, threshold = 0.5)
#'
#' @export
#'
#' @seealso \code{read_fasta_alignment}
remove_gaps <- function(x, threshold = 0.5) {
  if (!(threshold >= 0.5 && threshold < 1)) {
    stop(paste0(
      "'threshold' must be between 0.5 and <1. \n
      You've set it to ", threshold,"."
    ), call. = FALSE)
  }
  if (!is.rprimer_alignment(x)) {
    stop("'x' must be an rprimer_alignment object.", call. = FALSE)
  }
  # Catch sequence names
  sequence_names <- names(x)
  # Do a strsplit
  splitted <- purrr::map(x, split_sequence)
  # Make a matrix
  matr <- do.call("rbind", splitted)
  # gaps will be represented as 1
  matr <- gsub("-", 1, matr)
  # Other characters will be represented as 0
  matr <- gsub("[a-z]", 0, matr)
  # Convert to integer
  matr <- apply(matr, 2, as.integer)
  # Calculate gap frequency
  if (!is.matrix(matr)) {
    matr <- matrix(matr, ncol = length(matr)) # If only one sequence
  }
  gaps <- colMeans(matr)
  # Get the indexes of the positions that should be kept
  index <- which(gaps <= threshold)
  # Keep positions with gap frequency below or at threshold
  splitted <- purrr::map(seq_along(splitted), ~ splitted[[.x]][index])
  # And paste them together
  x <- purrr::map(splitted, paste, collapse = "")
  names(x) <- sequence_names
  x <- new_rprimer_alignment(x)
}
