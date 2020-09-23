#' Reverse complement
#'
#' \code{reverse_complement} finds the reverse complement of a DNA seuquence.
#'
#' @param x A DNA sequence (a character vector of length one).
#'
#' @details For \code{x}, valid bases are 'a', 'c', 'g', 't', 'r', 'y', 'm',
#' 'k', 's', 'w', n', 'h', 'd', 'v', 'b' and '-'.
#'
#' @return The reverse complement. Non valid bases will return as \code{NA}.
#'
#' @examples
#' reverse_complement("cttgtr")
#'
#' @export
reverse_complement <- function(x) {
  if (typeof(x) != "character") {
    stop("'x' must be a character vector", call. = FALSE)
  }
  x <- tolower(x)
  if (grepl("[^acgtrymkswnhdvb-]", x)) {
    stop("'x' contains at least one invalid base. \n
      Valid bases are 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w',
      'n', 'h', 'd', 'v', 'b' and '-'",
         call. = FALSE
    )
  }
  x <- strsplit(x, split = "")
  complement <- complement_lookup[unlist(x)]
  complement <- unname(complement)
  rc <- rev(complement)
  rc <- paste(rc, collapse = "")
  rc
}
