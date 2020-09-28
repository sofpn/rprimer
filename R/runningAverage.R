#' Calculate running average
#'
#' \code{runningAverage} calculates a running average of a
#' numeric vector.
#'
#' @param x A numeric vector.
#'
#' @param size
#' Optional. The number of observations in each average, which
#' can range from \code{1} to \code{length(x)}.
#' If not specified, the size will be set to the nearest positive, nonzero
#' integer to \code{length(x)/100}.
#'
#' @return A tibble with position and running average of x.
#'
#' @keywords internal
runningAverage <- function(x, size = NULL) {
    if (is.null(size)) {
        size <- round(length(x) / 100)
        if (size == 0) {
            size <- 1
        }
    }
    if (!is.numeric(size)) {
        stop("'size' must in numeric format.", call. = FALSE)
    }
    if (size > length(x)) {
        stop(
            paste0("'size' cannot be greater than ", length(x), "."),
            call. = FALSE
        )
    }
    if (size < 1) {
        stop("'size' must be greater than 1.", call. = FALSE)
    }
    sums <- c(0, cumsum(x))
    from <- seq_len(length(sums) - size)
    to <- seq(size + 1, length(sums))
    average <- (sums[to] - sums[from]) / size
    position <- seq(size, length(x))
    if (size > 1) {
        midpoint <- size / 2
        position <- position - midpoint
    }
    df <- tibble::tibble(position, average)
    df
}
