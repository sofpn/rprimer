#' Construct rprimer_profile objects
#'
#' @param x An rprimer_profile-like object.
#'
#' @details an rprimer_profile object must be a matrix of type
#' 'double'. Each element must have a value between 0 and 1. The matrix
#' must have both rownames and colnames.
#'
#' @return An rprimer_profile object if the validation is succeeds.
#' An error message if not.
#'
#' @keywords internal
#'
#' @noRd
new_rprimer_profile <- function(x = matrix()) {
  # Integrity checks
  stopifnot(
    is.matrix(x), is.double(x), max(x) <= 1, min(x) >= 0,
    !is.null(rownames(x)), !is.null(colnames(x))
  )
  # Set class attribute
  x <- structure(x, class = "rprimer_profile")
  return(x)
}

#' Check if an object is an rprimer_profile
#'
#' @param x An R object.
#'
#' @return \code{TRUE} or \code{FALSE}.
#'
#' @keywords internal
#'
#' @noRd
is.rprimer_profile <- function(x) {
  inherits(x, "rprimer_profile")
}
