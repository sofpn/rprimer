#' Get sequence properties
#'
#' \code{sequence_properties} returns sequence information from a sequence profile
#'
#' @param x A sequence profile (an object of class 'rprimer_profile').
#'
#' @param iupac_threshold
#' A number between greater than 0 and less or equal to 0.2.
#' At each position, all nucleotides with a proportion
#' higher than or equal to the stated threshold will be included in
#' the iupac consensus sequence. The default is \code{NULL},
#' which means that all nucleotides that are present at the
#' position in matter will be included.
#'
#' @section
#' Majority consensus sequence:
#' The most frequently occuring nucleotide.
#' If two or more bases occur with the same frequency,
#' the consensus nucleotide will be randomly selected among these bases.
#'
#' @section
#' IUPAC consensus sequence:
#' The consensus sequence expressed in iupac format (i.e. with wobble bases)
#' Note that the iupac consensus sequence only
#' takes 'a', 'c', 'g', 't' and '-' as input. Degenerate bases
#' present in the alignment will be skipped. If a position only contains
#' degenerate/invalid bases, the IUPAC consensus will be \code{NA} at that
#' position.
#'
#' @section Gaps:
#' Gaps are recognised as "-" in the sequence profile.
#'
#' @section Identity:
#' The nucleotide identity is the proportion of
#' the most common base. Gaps (-),
#' as well as nucleotides other than a, c, g and t are excluded from the
#' calculation.
#'
#' @section Entropy:
#' Shannon entropy is a measurement of
#' variability. First, for each nucleotide that occurs at a specific position,
#' \code{p*log2(p)}, is calculated, where \code{p} is the proportion of
#' that nucleotide. Then, the shannon entropy is calculated by summarising
#' these values for each nucleotide at the position in matter,
#' followed by multiplication by \code{-1}.
#' A value of \code{0} indicate no variability and a high value
#' indicate high variability.
#' Gaps (-), as well as bases other than
#' a, c, g and t are excluded from the calculation.
#'
#' @return
#' A tibble (data frame) of class 'rprimer_properties',
#' with information about majority and iupac consensus sequence, gap frequency,
#' nucleotide identity and shannon entropy.
#'
#' @examples
#' sequence_properties(example_rprimer_profile)
#' sequence_properties(example_rprimer_profile, iupac_threshold = 0.1)
#'
#' @export
sequence_properties <- function(x, iupac_threshold = NULL) {
  if (!is.rprimer_profile(x)) {
    stop("'x' must be an rprimer_profile object", call. = FALSE)
  }
  position <- seq_len(ncol(x))
  majority <- majority_consensus(x)
  iupac <- iupac_consensus(x, threshold = iupac_threshold)
  gaps <- gap_frequency(x)
  identity <- nucleotide_identity(x)
  entropy <- shannon_entropy(x)
  sequence_properties <- tibble::tibble(
    position, majority, iupac, gaps, identity, entropy
  )
  sequence_properties <- tibble::new_tibble(
    sequence_properties, nrow = nrow(sequence_properties),
    class = "rprimer_properties"
  )
  sequence_properties
}

#' Majority consensus sequence
#'
#' \code{majority_consensus} returns the majority consensus sequence of an
#' alignment of DNA sequences.
#'
#' @inheritParams sequence_properties
#'
#' @return The majority consensus sequence (a character vector of length n).
#'
#' @noRd
majority_consensus <- function(x) {
  if (!is.rprimer_profile(x)) {
    stop("'x' must be an rprimer_profile object.", call. = FALSE)
  }
  # Function to identify the most common base at a position
  find_most_common_base <- function(x, y) {
    most_common <- rownames(x)[which(y == max(y))]
    # If there are ties, the most common base will be randomly selected
    if (length(most_common > 1)) {
      most_common <- sample(most_common, 1)
    }
    return(most_common)
  }
  # Get the consensus sequence at all positions
  consensus <- apply(x, 2, function(y) find_most_common_base(x, y))
  consensus <- unname(consensus)
  consensus
}

#' Convert DNA nucleotides into the corresponding iupac degenerate base
#'
#' \code{as_iupac} takes several DNA nucleotides as input,
#' and returns the degenerate base in iupac format.
#'
#' @param x A character vector of length one, containing DNA
#' nucleotides (valid bases are a, c, g, t, or -). Each base must
#' be separated by a comma (,), e.g. 'a,c,g'.
#' Characters other than a, c, g, t and - will be ignored. However, when
#' the input only consist of invalid bases,
#' or if the bases are not separated by ',',
#' \code{as_iupac} will return NA.
#'
#' @return The corresponding iupac base.
#'
#' @examples
#' as_iupac("a,c")
#' as_iupac("r")
#' as_iupac("tg") # Will return NA since the bases are not separated by comma
#' @noRd
as_iupac <- function(x) {
  if (!(is.character(x) && length(x) == 1)) {
    stop(
      "'x' must be a character vector of length one, e.g. 'a,c,t'.",
      call. = FALSE
    )
  }
  x <- gsub(" ", "", x)
  x <- unlist(strsplit(x, split = ","), use.names = FALSE)
  x <- x[order(x)]
  x <- unique(x)
  bases <- c("a", "c", "g", "t", "-") # accepted nucleotides
  match <- x %in% bases
  x <- x[!(match == FALSE)] # exclude non accepted nucleotides
  x <- paste(x, collapse = ",")
  iupac <- unname(iupac_lookup[x]) # match with lookup table
  iupac
}

#' Iupac consensus sequence
#'
#' \code{iupac_consensus} returns the iupac consensus sequence from an
#' alignment of DNA sequences.
#'
#' @inheritParams sequence_properties
#'
#' @param threshold
#' A number between greater than 0 and less or equal to 0.2.
#' At each position, all nucleotides with a proportion
#' higher than or equal to the stated threshold will be included in
#' the iupac consensus sequence. The default is \code{NULL},
#' which means that all nucleotides that are present at the
#' position in matter will be included.
#'
#' @return The consensus sequence (a character vector of length n).
#'
#' @seealso \code{as_iupac} for further info on how the iupac consensus
#' sequecnce is determined.
#'
#' @noRd
iupac_consensus <- function(x, threshold = NULL) {
  if (!is.rprimer_profile(x)) {
    stop("'x' must be an rprimer_profile object.", call. = FALSE)
  }
  if (is.null(threshold)) {
    threshold <- 0
  }
  if (!is.double(threshold) || threshold < 0 || threshold > 0.2) {
    stop(paste0(
      "'treshold' must be between 0 and 0.2. \n
      You've set it to ", threshold, "."
    ))
  }
  # Select only a, c, g, t and - to count in the iupac consensus sequence
  bases <- c("a", "c", "g", "t", "-")
  x <- x[rownames(x) %in% bases, ]
  bases_to_include <- apply(x, 2, function(y) {
    paste(rownames(x)[y > threshold], collapse = ",")
  })
  bases_to_include <- unname(bases_to_include)
  consensus <- purrr::map_chr(bases_to_include, ~ as_iupac(.x))
  if (any(is.na(consensus))) {
    warning("The consensus sequence contain NAs. \n
    Try to lower the threshold value.", call. = FALSE)
  }
  consensus
}

#' Gap frequency
#'
#' \code{gap_frequency} returns the gap frequency at each position in an
#' alignment of DNA sequences.
#'
#' @inheritParams sequence_properties
#'
#' @return The gap frequency (a numeric vector of length n).
#'
#' @noRd
gap_frequency <- function(x) {
  if (!is.rprimer_profile(x)) {
    stop("'x' must be an rprimer_profile object.", call. = FALSE)
  }
  gaps <- x[rownames(x) == "-", ]
  gaps <- unname(gaps)
  gaps
}

#' Nucleotide identity
#'
#' \code{nucleotide_identity} returns the nucleotide identity from an
#' alignment of DNA sequences.
#'
#' @inheritParams sequence_properties
#'
#' @return The nucleotide identity (a numeric vector of length n).
#' The nucleotide identity can range from > 0 to 1.
#'
#' @noRd
nucleotide_identity <- function(x) {
  if (!is.rprimer_profile(x)) {
    stop("'x' must be an rprimer_profile object.", call. = FALSE)
  }
  # We want to assess identity based on DNA bases,
  # i.e. ignore gaps and degenerate positions, so we make a subset (s) of x
  # with the rows named a, c, g and t.
  bases <- c("a", "c", "g", "t")
  s <- x[rownames(x) %in% bases, ]
  # Calculate relative proportions of the bases in s
  s <- apply(s, 2, function(x) x / sum(x))
  # Find the largest proportion at each position
  identity <- apply(s, 2, max)
  identity <- unname(identity)
  identity[is.na(identity)] <- 0
  identity
}

#' Shannon entropy
#'
#' #' \code{shannon_entropy} returns the Shannon entropy from an
#' alignment of DNA sequences.
#'
#' @inheritParams sequence_properties
#'
#' @return The Shannon entropy (a numeric vector of length n).
#'
#' @noRd
shannon_entropy <- function(x) {
  if (!is.rprimer_profile(x)) {
    stop("'x' must be an rprimer_profile object.", call. = FALSE)
  }
  # We want to assess identity based on DNA bases,
  # i.e. ignore gaps and degenerate positions, so we make a subset (s) of x
  # with the rows named a, c, g and t.
  bases <- c("a", "c", "g", "t")
  s <- x[rownames(x) %in% bases, ]
  # Calculate relative proportions of the bases in s
  s <- apply(s, 2, function(x) x / sum(x))
  entropy <- apply(s[rownames(s) %in% bases, ], 2, function(x) {
    ifelse(x == 0, 0, x * log2(x))
  })
  entropy <- -colSums(entropy)
  entropy <- unname(entropy)
  entropy <- abs(entropy) # abs to avoid -0 (due to neg sums)
  entropy[is.na(entropy)] <- 0
  entropy
}

#' Check if an object has class attribute rprimer_properties
#'
#' @param x An 'rprimer_properties'-like object.
#'
#' @return \code{TRUE} or \code{FALSE}.
#'
#' @noRd
is.rprimer_properties <- function(x) inherits(x, "rprimer_properties")
