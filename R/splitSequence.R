#' Split sequence
#'
#' @param x A character vector of length one.
#'
#' @return A character vector of length \code{nchar(x)}.
#'
#' @example splitSequence("GCCTTG")
#'
#' @keywords internal
#'
#' @noRd
splitSequence <- function(x) {
    x <- unlist(strsplit(x, split = ""), use.names = FALSE)
    x
}
