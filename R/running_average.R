#' Calculate running average
#'
#' \code{running_average} calculates the centered running average of a
#' numeric vector.
#'
#' @param x A numeric vector.
#'
#' @param size The number of observations in each average.
#' If \code{NULL}, the size will be set to the nearest positive, nonzero
#' integer to \code{length(x)/100}.
#'
#' @return A tibble with position and running average of x.
#'
#' @examples
#' running_average(c(1, 3, 2, 1, 3, 4, 5), size = 2)
#' running_average(c(10, 22.4, 14.2, 44, 32))
#'
#' @keywords internal
#'
#' @noRd
running_average <- function(x, size = NULL) {
  if (is.null(size)) {
    size <- round(length(x) / 100)
    if (size == 0) {
      size <- 1
    }
  }
  stopifnot(is.numeric(x) && size <= length(x) && size >= 1)
  sums <- c(0, cumsum(x))
  average <- (sums[((size + 1):length(sums))] - sums[(1:(length(sums) - size))]) / size
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
