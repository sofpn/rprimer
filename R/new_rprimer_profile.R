# When implementing a matrix/array class, you should implement these methods:
# dim (gets you nrow and ncol), t, dimnames (gets you rownames and colnames),
# dimnames<- (gets you colnames<-, rownames<-), cbind, rbind.


#' Construct, subset and check 'rprimer_profile' objects
#'
#' * 'new_rprimer_profile()' constructs and validate an rprimer_profile-
#'     object
#' * 'is_rprimer_profile()' checks if an object has class attribute
#'     'rprimer_profile'
#' * '`[.rprimer_profile()`' subsets an object of class 'rprimer_profile'
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

#' @describeIn new_rprimer_profile
#'
#' @param x An rprimer_profile object.
#'
#' @return \code{TRUE} or \code{FALSE}.
#'
#' @noRd
is.rprimer_profile <- function(x) inherits(x, "rprimer_profile")

