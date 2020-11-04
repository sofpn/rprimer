#' Calculate running average
#'
#' \code{.runningAverage} calculates a running average of a
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
#' @return A data frame with position and running average of x.
#'
#' @keywords internal
#'
#' @noRd
.runningAverage <- function(x, size = NULL) {
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
    data.frame(position, average)
}

#' Calculate GC running average
#'
#' \code{.gcRunningAverage} calculates a running average of
#' the GC-content of a DNA sequence
#'
#' @param x
#' A DNA sequence (a character vector). Valid bases are ACGT-.
#'
#' @param size
#' Optional. The number of observations in each average, which
#' can range from \code{1} to \code{length(x)}.
#' If not specified, the size will be set to the nearest positive, nonzero
#' integer to \code{length(x)/100}.
#'
#' @return A data frame with position and running average of x.
#'
#' @seealso gcContent
#'
#' @keywords internal
#'
#' @noRd
.gcRunningAverage <- function(x, size = NULL) {
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
    x <- toupper(x)
    begins <- seq_len(length(x) - size + 1)
    ends <- begins + size - 1
    average <- purrr::map2_dbl(begins, ends, function(i, j) {
        string <- paste(x[i:j], collapse = "")
        .gcContent(string)
    })
    position <- seq(size, length(x))
    if (size > 1) {
        midpoint <- size / 2
        position <- position - midpoint
    }
    data.frame(position, average)
}
