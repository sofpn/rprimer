#' Calculate GC content of a DNA sequence
#'
#' \code{gcContent} finds the GC content of a DNA sequence.
#'
#' @param x
#' A DNA sequence (a character vector of length one).
#' Valid bases are ACGT-.
#'
#' @return The GC content of x.
#'
#' @keywords internal
gcContent <- function(x) {
  if (typeof(x) != "character" || length(x) != 1) {
    stop("'x' must be a character vector of length one.", call. = FALSE)
  }
  x <- toupper(x)
  if (grepl(paste0("[^", dnaBases, "]"), x)) {
    stop(paste0("'x' can only contain bases ", dnaBases, "."),
         call. = FALSE
    )
  }
  x <- splitSequence(x)
  gcCount <- length(which(x == "C" | x == "G"))
  totalCount <- length(which(x == "A" | x == "C" | x == "G" | x == "T"))
  gc <- gcCount / totalCount
  gc
}
