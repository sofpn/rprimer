#' Calculate running average of GC content
#'
#' \code{gc_running_average} calculates the centered running average of
#' the GC-content of a DNA sequence
#'
#' @param x a DNA sequence (a character vector of length one).
#'
#' @param size The number of observations in each average.
#' If \code{NULL}, the size will be set to the nearest positive, nonzero
#' integer to \code{length(x)/100}.
#'
#' @details \code{x} cannot contain other characters than
#' 'a', 'c', 'g', 't' and '-'.
#'
#' @return A tibble with position and running average of x.
#'
#' @examples
#' running_average(c("gtttgtttctt"), size = 2)
#' running_average(c("cttggggtttctttggtt-ttagg"))
#'
#' @seealso
#' gc_content
#'
#' @keywords internal
#'
#' @noRd
gc_running_average <- function(x, size = NULL) {
  if (is.null(size)) {
    size <- round(length(x) / 100)
    if (size == 0) {
      size <- 1
    }
  }
  stopifnot(is.character(x) && size <= length(x) && size >= 1)
  begins <- seq_len(length(x) - size + 1)
  ends <- begins + size - 1
  # Calculate GC content of chunks of x
  average <- purrr::map2_dbl(begins, ends, function(y, z) {
    string <- paste(x[y:z], collapse = "")
    gc_content(string)
  })
  # Set the position to get a trailed moving average
  position <- seq(size, length(x))
  # To get a centered running average,
  # I align each moving average to the midpoint of the range of observations
  if (size > 1) {
    midpoint <- size / 2
    position <- position - midpoint
  }
  df <- tibble::tibble(position, average)
  df
}
