#' Construct, subset and check 'rprimer_alignment' objects
#'
#' * 'new_rprimer_alignment()' constructs and validate an rprimer_alignment-
#'     object
#' * 'is_rprimer_alignment()' checks if an object has class attribute
#'     'rprimer_alignment'
#' * '`[.rprimer_alignment()`' subsets an object of class 'rprimer_alignment'
#'
#' @param x An rprimer_alignment-like object.
#'
#' @details
#' An rprimer_alignment object must must be a list and
#' contain at least one DNA
#' sequence, and all sequences must be a character vector of lentgh one,
#' in lowercase format.
#' All sequences (including gaps) must be of the same length.
#' All sequences must have unique names. Valid bases are
#' a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w',
#' 'n', 'h', 'd', 'v', 'b' and '-')
#'
#' @return
#' An rprimer_alignment object if the validation is succeeds.
#' An error message if not.
#'
#' @noRd
new_rprimer_alignment <- function(x = list()) {
  # x must be a list
  stopifnot(is.list(x))
  # All sequences must be a character vector of length one
  has_length_one <- purrr::map_lgl(x, ~ is.character(.x) && length(.x) == 1)
  if (any(has_length_one == FALSE)) {
    stop(
      "All sequences must be a character vector of length one. \n
      At least one sequence is not a character vector/is not of length one.",
      call. = FALSE
    )
  }
  # All sequences must contain the same number of characters
  sequence_lengths <- purrr::map_int(x, nchar)
  if (length(unique(sequence_lengths)) != 1) {
    stop(
      "The sequences does not appear to be aligned. \n
      All sequences (including gaps) are not of the same length.",
      call. = FALSE
    )
  }
  # All sequence names must be unique
  unique_name <- length(unique(names(x))) == length(x)
  if (unique_name == FALSE) {
    stop("All sequences must have unique names.", call. = FALSE)
  }
  # All sequences must be in lowercase format and contain only valid bases
  non_valid_base <- purrr::map_lgl(x, ~ grepl("[^acgtrymkswnhdvb-]", .x))
  if (any(non_valid_base)) {
    stop(
      "At least one sequence contain one or more invalid characters. \n
      Valid characters are: \n
      'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w', 'n', 'h', 'd',
      'v', 'b' and '-'",
      call. = FALSE
    )
  }
  # Set class attribute
  x <- structure(x, class = "rprimer_alignment")
  x
}

#' @describeIn new_rprimer_alignment
#'
#' @param x An rprimer_alignment object.
#'
#' @retun \code{TRUE} or \code{FALSE}.
#'
#' @noRd
is.rprimer_alignment <- function(x) inherits(x, "rprimer_alignment")

#' @describeIn new_rprimer_alignment
#'
#' @param x An rprimer_alignment object.
#'
#' @return A subset.
#'
#' @noRd
`[.rprimer_alignment` <- function(x, i, ...) {
  new_rprimer_alignment(NextMethod())
}
