#' Split sequence
#'
#' @param x A character vector of length one.
#'
#' @return A character vector of length \code{nchar(x)}.
#'
#' @example split_sequence("cgtttg")
#'
#' @noRd
split_sequence <- function(x) {
  stopifnot(is.character(x), length(x) == 1)
  x <- unlist(strsplit(x, split = ""), use.names = FALSE)
  x
}

#' Truncate the name of a sequence in fasta format
#'
#' \code{truncate_name} shortens the name of a sequence.
#'
#' @param x A sequence name (a character vector of length one).
#'
#' @return A shorter version of the sequence name (the first
#' word of the sequence name). '>' symbols are removed.
#'
#' @example
#' truncate_name('AB856243.1 Hepatitis E virus')
#'
#' @noRd
truncate_name <- function(name) {
  name <- strsplit(name, split = "[[:space:]]")
  name <- unlist(name, use.names = FALSE)
  name <- name[[1]]
  name <- gsub(">", "", name)
  name
}

#' Complement
#'
#' \code{complement} finds the complement of a DNA sequence.
#'
#' @param x a DNA sequence (a character vector of length one, e.g. 'cttgg').
#'
#' @details For \code{x}, valid bases are 'a', 'c', 'g', 't', 'r', 'y', 'm',
#' 'k', 's', 'w', n', 'h', 'd', 'v', 'b' and '-'.
#'
#' @return The complement sequence of x.
#'
#' @examples
#' reverse_complement("cttgtr")
#' @noRd
complement <- function(x) {
  if (typeof(x) != "character") {
    stop("'x' must be a character vector.", call. = FALSE)
  }
  x <- tolower(x)
  if (grepl("[^acgtrymkswnhdvb-]", x)) {
    stop("'x' contains at least one invalid base.
      Valid bases are 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w',
      'n', 'h', 'd', 'v', 'b' and '-'",
      call. = FALSE
    )
  }
  x <- strsplit(x, split = "")
  complement <- complement_lookup[unlist(x)]
  complement <- unname(complement)
  complement
}

