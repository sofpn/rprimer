#' Reverse complement
#'
#' \code{reverseComplement} finds the reverse complement of a DNA sequence.
#'
#' @param x
#' A DNA sequence (a character vector of length one). Valid bases are
#' ACGTRYSWKMBDHVN-.
#'
#' @return
#' The reverse complement of x. Non valid bases will return as \code{NA}.
#'
#' @examples
#' reverseComplement("CTTTGRTN")
#' reverseComplement("CTTTGRTN-")
#' reverseComplement("C")
#'
#' @keywords internal
reverseComplement <- function(x) {
  if (typeof(x) != "character" || length(x) != 1) {
    stop("'x' must be a character vector of length one.", call. = FALSE)
  }
  x <- toupper(x)
  if (grepl(paste0("[^", allBases, "]"), x)) {
    stop(paste0("'x' can only contain bases ", allBases, "."),
         call. = FALSE
    )
  }
  x <- splitSequence(x)
  complement <- complementLookup[unlist(x)]
  complement <- unname(complement)
  rc <- rev(complement)
  rc <- paste(rc, collapse = "")
  rc
}
