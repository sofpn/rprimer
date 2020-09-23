#' Calculate GC running average
#'
#' \code{gcRunningAverage} calculates running average of
#' the GC-content of a DNA sequence
#'
#' @param x a DNA sequence (a character vector).
#'
#' @param size Optional. The number of observations in each average.
#' can range from \code{1} to \code{length(x)}.
#' If not specified, the size will be set to the nearest positive, nonzero
#' integer to \code{length(x)/100}.
#'
#' @details \code{x} cannot contain other characters than
#' 'A', 'C', 'G', 'T' and '-'.
#'
#' @return A tibble with position and running average of x.
#'
#' @examples
#' gcRunningAverage(
#' c("G", "T", "T", "T", "G", "T", "T", "T", "C", "T", "T"), size = 2
#' )
#' @seealso
#' gcContent
#'
#' @keywords internal
#'
#' @noRd
gcRunningAverage <- function(x, size = NULL) {
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
  if (grepl(paste0("[^", dnaBases, "]"), paste(x, collapse = ""))) {
    stop(paste0("'x' can only contain bases ", dnaBases, "."),
         call. = FALSE
    )
  }
  begins <- seq_len(length(x) - size + 1)
  ends <- begins + size - 1
  average <- purrr::map2_dbl(begins, ends, function(i, j) {
    string <- paste(x[i:j], collapse = "")
    gcContent(string)
  })
  position <- seq(size, length(x))
  if (size > 1) {
    midpoint <- size / 2
    position <- position - midpoint
  }
  df <- tibble::tibble(position, average)
  df
}
