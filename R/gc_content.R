#' Calculate GC content of a DNA sequence
#'
#' \code{gc_content} finds the GC content of a DNA sequence.
#'
#' @param x a DNA sequence (a character vector of length one).
#'
#' @details \code{x} cannot contain other characters than
#' 'a', 'c', 'g', 't' and '-'.
#'
#' @return the GC content of x. Gaps ('-') will not be included
#' in the calculation.
#'
#' @examples
#' gc_content("acgttcc")
#' gc_content("acgttcc--")
#' gc_content("acgrn") ## Will return an error because of an invalid base.
#'
#' @export
gc_content <- function(x) {
  if (typeof(x) != "character" || length(x) != 1) {
    stop("'x' must be a character vector of length one.", call. = FALSE)
  }
  x <- tolower(x)
  if (grepl("[^acgt-]", x)) {
    stop("'x' contains at least one invalid base. \n
      Valid bases are 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w',
      'n', 'h', 'd', 'v', 'b' and '-'",
         call. = FALSE
    )
  }
  x <- split_sequence(x)
  gc_count <- length(which(x == "c" | x == "g"))
  # Gaps will not be included in the total count
  total_count <- length(which(x == "a" | x == "c" | x == "g" | x == "t"))
  gc <- gc_count / total_count
  gc
}
