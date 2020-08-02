# Utils =======================================================================

#' Split sequence
#'
#' @param x A character vector of length one
#'
#' @return A character vector of length \code{nchar(x)}
#'
#' @example split_sequence("cgtttg")
#'
#' @noRd
split_sequence <- function(x) {
  stopifnot(is.character(x), length(x) == 1)
  x <- unlist(strsplit(x, split = ""), use.names = FALSE)
  return(x)
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
  return(name)
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
  if (typeof(x) != "character" || length(x) != 1) {
    stop("x must be a character vector of length one", call. = FALSE)
  }
  x <- tolower(x)
  if (grepl("[^acgtrymkswnhdvb-]", x)) {
    stop("x contains at least one invalid base.
      Valid bases are 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w',
      'n', 'h', 'd', 'v', 'b' and '-'",
      call. = FALSE
    )
  }
  x <- strsplit(x, split = "")
  complement <- complement_lookup[unlist(x)]
  complement <- unname(complement)
  return(complement)
}

# Context: Get sequence properties ============================================

#' Majority consensus sequence
#'
#' \code{majority_consensus} returns the majority consensus sequence of an
#' alignment of DNA sequences.
#'
#' @param x A nucleotide profile of class 'rprimer_profile', i.e.
#' a numeric m x n matrix, where m is the number of nucleotides in the
#' alignment, and n is the number of positions in the alignment.
#'
#' @details If there are ties (two or more bases occur with the same frequency)
#' , the consensus base will be randomly selected among these bases.
#'
#' @return the majority consensus sequence (a character vector of length n).
#'
#' @noRd
majority_consensus <- function(x) {
  if (!inherits(x, "rprimer_profile")) {
    stop("An rprimer_profile object is expected.", call. = FALSE)
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
  return(consensus)
}

#' Convert DNA nucleotides into the corresponding iupac degenerate base
#'
#' \code{as_iupac} takes several DNA nucleotides as input,
#' and returns the degenerate base in iupac format.
#'
#' @param x a character vector of length one containing DNA
#' nucleotides (valid bases are a, c, g, t, or -). Each base must
#' be separated by a comma (,), e.g. 'a,c,g'.
#' Characters other than a, c, g, t and - will be ignored. However, when
#' the input only consist of invalid bases,
#' or if the bases are not separated by ',',
#' \code{as_iupac} will return NA.
#'
#' @return the corresponding iupac base (a character vector of length one).
#'
#' @examples
#' as_iupac("a,c")
#' as_iupac("r")
#' as_iupac("tg") # Will return NA since the bases are not separated by comma
#' @noRd
as_iupac <- function(x) {
  if (!(is.character(x) && length(x) == 1)) {
    stop(
      "A character vector of length one is expected, e.g. 'a,c,t'",
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
  return(iupac)
}

#' Iupac consensus sequence
#'
#' \code{iupac_consensus} returns the iupac consensus sequence from an
#' alignment of DNA sequences.
#'
#' @param x A nucleotide profile of class 'rprimer_profile', i.e.
#' a numeric m x n matrix, where m is the number of nucleotides in the
#' alignment, and n is the number of positions in the alignment.
#'
#' @param threshold A number between 0 and 0.2 (the default is 0).
#' At each position, all nucleotides with a proportion higher than
#' or equal to the stated threshold will be included in the iupac consensus
#' sequence.
#'
#' @details Note that \code{iupac_consensus} only
#' takes 'a', 'c', 'g', 't' and '-' as input. Degenerate bases
#' present in the alignment will be skipped. If a position only contains
#' degenerate/invalid bases the iupac consensus character will be NA at that
#' position.
#'
#' @return The consensus sequence (a character vector of length n).
#'
#' @seealso \code{as_iupac} for further info on how the iupac consensus
#' sequecnce is determined.
#'
#' @noRd
iupac_consensus <- function(x, threshold = 0) {
  if (!inherits(x, "rprimer_profile")) {
    stop("An rprimer_profile object is expected.", call. = FALSE)
  }
  if (!is.double(threshold) || threshold < 0 || threshold > 0.2) {
    stop(paste0(
      "The threshold was set to ", threshold, ". Valid threshold values
     are from 0 to 0.2"
    ))
  }
  # Select only a, c, g, t and - to count in the iupac consensus sequence
  bases <- c("a", "c", "g", "t", "-")
  x <- x[rownames(x) %in% bases, ]
  bases_to_include <- apply(x, 2, function(y) {
    paste(rownames(x)[y >= threshold], collapse = ",")
  })
  bases_to_include <- unname(bases_to_include)
  consensus <- purrr::map_chr(bases_to_include, ~ as_iupac(.x))
  if (any(is.na(consensus))) {
    warning("The consensus sequence contain NAs.
    Try to lower the threshold value.", call. = FALSE) # Check if this works!
  }
  return(consensus)
}

#' Gap frequency
#'
#' \code{gap_frequency} returns the gap frequency at each position in an
#' alignment of DNA sequences.
#'
#' @param x A nucleotide profile of class 'rprimer_profile', i.e.
#' a numeric m x n matrix, where m is the number of nucleotides in the
#' alignment, and n is the number of positions in the alignment.
#'
#' @details Gaps are recognised as "-".
#'
#' @return the gap frequency (a numeric vector of length n).
#'
#' @noRd
gap_frequency <- function(x) {
  if (!inherits(x, "rprimer_profile")) {
    stop("An rprimer_profile object is expected.", call. = FALSE)
  }
  gaps <- x[rownames(x) == "-", ]
  gaps <- unname(gaps)
  return(gaps)
}

#' Nucleotide identity
#'
#' #' \code{nucleotide_identity} returns the nucleotide identity from an
#' alignment of DNA sequences.
#'
#' @param x A nucleotide profile of class 'rprimer_profile', i.e.
#' a numeric m x n matrix, where m is the number of nucleotides in the
#' alignment, and n is the number of positions in the alignment.
#'
#' @details The nucleotide identity is the proportion of
#' the most common base at each position in the alignment.  Gaps (-),
#' as well as nucleotides other than a, c, g and t, are excluded from the
#' calculation.
#'
#' @return the nucleotide identity (a numeric vector of length n).
#' The nucleotide identity can range from > 0 to 1.
#'
#' @noRd
nucleotide_identity <- function(x) {
  if (!inherits(x, "rprimer_profile")) {
    stop("An rprimer_profile object is expected.", call. = FALSE)
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
  return(identity)
}

#' Shannon entropy
#'
#' #' \code{shannon_entropy} returns the shannon entropy from an
#' alignment of DNA sequences.
#'
#' @param x A nucleotide profile of class 'rprimer_profile', i.e.
#' a numeric m x n matrix, where m is the number of nucleotides in the
#' alignment, and n is the number of positions in the alignment.
#'
#' @details Shannon entropy is a measurement of
#' variability, and it is here calculated at each position in the alignment.
#' First, for each nucleotide that occurs at the position in matter,
#'  \code{p*log2(p)}, is calculated, where p is the proportion of
#' that nucleotide. Then, the shannon entropy is calculated by summarising
#' these values for each nucleotide at the position in matter,
#' followed by multiplication by -1.
#' A value of 0 indicate no variability and a high value
#' indicate high variability.
#' Gaps (-), as well as nucleotides other than
#' a, c, g and t, are excluded from the calculation.
#'
#' @return the shannon entropy (a numeric vector of length n).
#'
#' @noRd
shannon_entropy <- function(x) {
  if (!inherits(x, "rprimer_profile")) {
    stop("An rprimer_profile object is expected.", call. = FALSE)
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
  return(entropy)
}

# Context: Get oligos =========================================================

#' Divide a DNA sequence into n-sized chunks
#'
#' \code{get_nmers} divides a character vector into chunks of size \code{n},
#' in steps of one.
#'
#' @param x A character vector.
#'
#' @param n The desired size of each 'chunk'/'mer'. An integer between
#' \code{1} and \code{length(x)}. The default is \code{NULL}.
#' In that case, n will be set to
#' the nearest integer to \code{length(x)/10}. However, if the nearest integer
#' is zero, \code{n} will be set to one.
#'
#' @return A character vector where each element is a 'mer' of size n.
#'
#' @examples
#' get_nmers(c("c", "g", "t", "t", "c", "g"), n = 2)
#' @noRd
get_nmers <- function(x, n = NULL) {
  if (!is.character(x)) {
    stop("x must be a character vector.", call. = FALSE)
  }
  if (is.null(n)) {
    n <- round(length(x) / 10)
    if (n == 0) {
      n <- 1
    }
  }
  if (!is.numeric(n) || n < 1 || n > length(x)) {
    stop("n must be an integer between 1 and length(x)", call. = FALSE)
  }
  begin <- 1:(length(x) - n + 1)
  end <- begin + n - 1
  nmer <- purrr::map_chr(
    begin, ~ paste(x[begin[[.x]]:end[[.x]]], collapse = "")
  )
  return(nmer)
}

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
#' gc_content("acgrn") # Will return an error because of an invalid base.
#' @noRd
gc_content <- function(x) {
  if (typeof(x) != "character" || length(x) != 1) {
    stop("x must be a character vector of length one", call. = FALSE)
  }
  x <- tolower(x)
  if (grepl("[^acgt-]", x)) {
    stop("x contains at least one invalid base.
      x can only contain 'a', 'c', 'g', 't' and '-'", call. = FALSE)
  }
  x <- split_sequence(x)
  gc_count <- length(which(x == "c" | x == "g"))
  # Gaps will not be included in the total count
  total_count <- length(which(x == "a" | x == "c" | x == "g" | x == "t"))
  gc <- gc_count / total_count
  return(gc)
}

#' Reverse complement
#'
#' \code{reverse_complement} finds the reverse complement of a DNA seuquence.
#'
#' @param x a DNA sequence (a character vector of length one, e.g. 'cttgg').
#'
#' @details For \code{x}, valid bases are 'a', 'c', 'g', 't', 'r', 'y', 'm',
#' 'k', 's', 'w', n', 'h', 'd', 'v', 'b' and '-'.
#'
#' @return The reverse complement. Non valid bases will return as NA.
#'
#' @examples
#' reverse_complement("cttgtr")
#' @noRd
reverse_complement <- function(x) {
  if (typeof(x) != "character" || length(x) != 1) {
    stop("x must be a character vector of length one", call. = FALSE)
  }
  x <- tolower(x)
  if (grepl("[^acgtrymkswnhdvb-]", x)) {
    stop("x contains at least one invalid base.
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
  return(rc)
}

#' Calculate running, cumulative sums
#'
#' \code{running_sum} calculates 'running' sums within a numeric vector. Each
#' sum is calculated in a size of \code{n}, in steps of 1 (i.e., if
#' \code{n = 20}, the sum will be calculated from element 1 to 20,
#' then from element 2 to 21, then from element 3 to 22, etc.)
#'
#' @param x A numeric vector.
#'
#' @param n The size of each sum that is to be calculated, an integer between
#' 1 and \code{length(x)}.
#'
#' @details The default of \conde{n} is \code{NULL}. In that case,
#' \code{n} will be set to
#' the nearest integer to \code{length(x)/10}. However, if
#' \code{length(x)/10} is zero, \code{n} will be
#' set to 1.
#'
#' @return The running sums of \code{x}, in steps of 1, in size of \code{n}
#' (a numeric vector of length \code{length(x) - n + 1}).
#'
#' @examples
#' running_sum(runif(100))
#' @noRd
running_sum <- function(x, n = NULL) {
  if (!(is.numeric(x))) {
    stop("A numeric vector is expected for x.", call. = FALSE)
  }
  if (is.null(n)) {
    n <- round(length(x) / 10)
    if (n == 0) {
      n <- 1
    }
  }
  if (!is.numeric(n) || n < 1 || n > length(x)) {
    stop("n must be a number between 1 and length(x)", call. = FALSE)
  }
  cumul <- c(0, cumsum(x))
  runsum <- cumul[(n + 1):length(cumul)] - cumul[1:(length(cumul) - n)]
  return(runsum)
}

#' Exclude oligos with a specific pattern
#'
#' \code{exclude} replaces oligos with a specific pattern
#' with \code{NA}.
#'
#' @param x One or more oligo sequences (a character vector).
#'
#' @param pattern A regular expression.
#'
#' @return A character vector where the oligos with the pattern
#' have been replaced with \code{NA}.
#'
#' @examples
#' exclude(c("cttgttatttt", "cgattctg"), "a$|t$")
#' @noRd
exclude <- function(x, pattern) {
  if (!is.character(x)) {
    stop("x must be a character vector", call. = FALSE)
  }
  if (!is.character(pattern) || length(pattern) != 1) {
    stop("pattern must be a character vector of length one", call. = FALSE)
  }
  regex <- pattern
  x[grepl(regex, x)] <- NA
  return(x)
}

#' Exclude oligos
#'
#' \code{exclude oligos} replaces oligos with many consecutive
#' mono- or dinucleotides with \code{NA}. If selected, it also replaces oligos
#' with 5'-end g or 3'-end t or a with \code{NA}.
#'
#' @param x One or more DNA sequences (a character vector).
#'
#' @param avoid_3end_ta \code{TRUE} or \code{FALSE}.
#' If oligos with t or a at the 3' end
#' should be replaced with \code{NA}. The default is \code{FALSE}.
#'
#' @param avoid_5end_g \code{TRUE} or \code{FALSE}.If oligos with g
#' at the 5' end should be replaced with \code{NA}. The default is \code{FALSE}.
#'
#' @param avoid_3end_runs \code{TRUE} or \code{FALSE}.
#' If oligos with more than two runs
#' of the same nucleotide at the 3' end should be replaced with \code{NA}.
#' The default is \code{FALSE}.
#'
#' @details An oligo is replaced with \code{NA} if
#' - It has more than three runs of the same dinucleotide (e.g. 'tatatata')
#' - It has more than four runs of the same nucleotide
#'
#' @return A character vector where unwanted oligos
#' has been replaced with \code{NA}.
#'
#' @examples
#' exclude_oligos(c("gtgtaaaaatatttt", "cgattctg"))
#' @noRd
exclude_oligos <- function(
                           x, avoid_3end_ta = FALSE, avoid_5end_g = FALSE, avoid_3end_runs = FALSE) {
  if (any(!is.logical(c(avoid_3end_ta, avoid_5end_g, avoid_3end_runs)))) {
    stop(
      "avoid_3end_ta, avoid_5end_g and avoid_3end_runs
        must be set to TRUE or FALSE",
      call. = FALSE
    )
  }
  # Remove oligos with at least 4 'runs' of the same dinucleotide
  x <- exclude(x, "(at|ta|ac|ca|ag|ga|gt|tg|cg|gc|tc|ct)\\1\\1\\1")
  # Remove oligos with at least 5 'runs' of the same nucleotide
  x <- exclude(x, "([a-z])\\1\\1\\1\\1")
  if (avoid_3end_ta == TRUE) {
    # Remove oligos with t or a at the 3' end
    x <- exclude(x, "a$|t$")
  }
  if (avoid_5end_g == TRUE) {
    # Remove oligos with g at the 5' end
    x <- exclude(x, "^g")
  }
  if (avoid_3end_runs == TRUE) {
    # Remove oligos with at least 3 'runs' of the same nucleotide in 3'
    x <- exclude(x, "([a-z])\\1\\1$")
  }
  return(x)
}

#' Exclude unwanted oligos
#'
#' \code{exclude unwanted oligos} replaces oligos with many consecutive
#' mono- or dinucleotides with \code{NA}. If selected, it also replaces oligos
#' with 5'-end g or 3'-end t or a with \code{NA}.
#'
#' @param x A tibble with oligos.
#'
#' @param avoid_3end_ta \code{TRUE} or \code{FALSE}.
#' If oligos with t or a at the 3' end
#' should be replaced with \code{NA}.
#'
#' @param avoid_5end_g \code{TRUE} or \code{FALSE}.If oligos with g
#' at the 5' end should be replaced with \code{NA}.
#'
#' @param avoid_3end_runs \code{TRUE} or \code{FALSE}.
#' If oligos with more than two runs
#' of the same nucleotide at the 3' end should be replaced with \code{NA}.
#'
#' @details An oligo is replaced with \code{NA} if
#' - It has more than three runs of the same dinucleotide (e.g. 'tatatata')
#' - It has more than four runs of the same nucleotide
#'
#' @return A tibble where unwanted oligos
#' has been replaced with \code{NA}.
#'
#' @noRd
exclude_unwanted_oligos <- function(
                                    x, avoid_3end_ta, avoid_5end_g, avoid_3end_runs) {
  x$majority <- exclude_oligos(
    x$majority,
    avoid_3end_ta = avoid_3end_ta,
    avoid_5end_g = avoid_5end_g,
    avoid_3end_runs = avoid_3end_runs
  )
  x$majority_rc <- exclude_oligos(
    x$majority_rc,
    avoid_3end_ta = avoid_3end_ta,
    avoid_5end_g = avoid_5end_g,
    avoid_3end_runs = avoid_3end_runs
  )
  x$iupac[is.na(x$majority)] <- NA
  x$iupac_rc[is.na(x$majority_rc)] <- NA
  # Identify oligos where both the sense and antisense sequence is NA
  invalid_oligos <- is.na(x$majority) & is.na(x$majority_rc)
  # Remove them
  x <- x[!invalid_oligos, ]
  return(x)
}

#' Count the number of degenerate bases in a DNA sequence
#'
#' \code{count_degenerates} returns the number of degenerate bases in a DNA
#' sequence.
#'
#' @param x a DNA sequence (a character vector of length one, e.g. 'cttgg').
#'
#' @details Valid bases for \code{x} are 'a', 'c', 'g', 't', 'r', 'y', 'm',
#' 'k', 's', 'w', n', 'h', 'd', 'v', 'b' and '-'.
#'
#' @return the number of degenerate bases in \code{x} (an integer).
#'
#' @examples
#' count_degenerates("cttnra")
#' @noRd
count_degenerates <- function(x) {
  if (typeof(x) != "character" || length(x) != 1) {
    stop("x must be a character vector of length one", call. = FALSE)
  }
  if (grepl("[^acgtrymkswnhdvb-]", x)) {
    stop("x contains at least one invalid base.
      Valid bases are 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w',
      'n', 'h', 'd', 'v', 'b' and '-'",
      call. = FALSE
    )
  }
  nt <- c("a", "c", "g", "t", "-")
  x <- split_sequence(x)
  count <- length(x[!x %in% nt])
  return(count)
}

#' Count the degeneracy of a DNA sequence
#'
#' \code{count_degenerates} counts the number of unique sequences of
#' a DNA sequence with degenerate bases.
#'
#' @param x a DNA sequence (a character vector of length one, e.g. 'cttgg').
#'
#' @details Valid bases for \code{x} are 'a', 'c', 'g', 't', 'r', 'y', 'm',
#' 'k', 's', 'w', n', 'h', 'd', 'v', 'b' and '-'.
#'
#' @return the number of unique sequences of x (an integer).
#'
#' @examples
#' count_degeneracy("cttnra")
#' @noRd
count_degeneracy <- function(x) {
  if (typeof(x) != "character" || length(x) != 1) {
    stop("x must be a character vector of length one", call. = FALSE)
  }
  if (grepl("[^acgtrymkswnhdvb-]", x)) {
    stop("x contains at least one invalid base.
      Valid bases are 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w',
      'n', 'h', 'd', 'v', 'b' and '-'",
      call. = FALSE
    )
  }
  x <- split_sequence(x)
  # Find the number of nucleotides at each position in x
  n_nucleotides <- degeneracy_lookup[x]
  # Calculate the total number of DNA sequences in x
  degeneracy <- prod(n_nucleotides)
  return(degeneracy)
}

#' Generate oligos of a specific length
#'
#' @param x An object of class 'rprimer_properties'
#'
#' @param oligo_length An integer. The minimum allowed
#' value is 6 and the maximum allowed value is 30. The default is 20.
#'
#' @param max_gap_frequency Maximum allowed gap frequency.
#' A number between 0 and 1 (default is 0.1, which means that only
#' positions with a gap frequency equal to or less than 0.1 will be
#' considered as an oligo region).
#'
#' @return A tibble with all possible oligos
#'
#' @noRd
generate_oligos <- function(
                            x,
                            oligo_length = 20,
                            max_gap_frequency = 0.1,
                            max_degenerates = 2,
                            max_degeneracy = 4) {
  if (!inherits(x, "rprimer_properties")) {
    stop(
      "An rprimer_properties object is expected for x.",
      call. = FALSE
    )
  }
  if (!(min(oligo_length) >= 6 && max(oligo_length) <= 30)) {
    stop("oligo_length must be between 6 and 30", call. = FALSE)
  }
  if (!(max_gap_frequency >= 0 && max_gap_frequency <= 1)) {
    stop("max_gap_frequency must be between 0 and 1", call. = FALSE)
  }
  if (!(max_degenerates <= 6 && max_degenerates >= 0)) {
    stop("max_degenerates must be between 0 and 6", call. = FALSE)
  }
  if (!(max_degeneracy >= 1 && max_degeneracy <= 64)) {
    stop("max_degeneracy must be between 1 and 64", call. = FALSE)
  }
  # Find all possible oligos of length y
  majority <- get_nmers(x$majority, n = oligo_length)
  iupac <- get_nmers(x$iupac, n = oligo_length)
  majority_rc <- purrr::map_chr(majority, ~ reverse_complement(.x))
  iupac_rc <- purrr::map_chr(iupac, ~ reverse_complement(.x))
  degenerates <- purrr::map_int(iupac, ~ count_degenerates(.x))
  degeneracy <- purrr::map_dbl(iupac, ~ count_degeneracy(.x))
  begin <- seq_along(majority)
  end <- seq_along(majority) + oligo_length - 1
  length <- oligo_length

  # Identify oligos with high gap frequency
  gap_bin <- ifelse(x$gaps > max_gap_frequency, 1L, 0L)
  gap_penalty <- running_sum(gap_bin, n = oligo_length)

  oligos <- tibble::tibble(
    begin, end, length, majority, iupac,
    majority_rc, iupac_rc, degenerates, degeneracy
  )
  # Exclude oligos with too high gap frequency
  oligos <- oligos[gap_penalty == 0, ]
  # Exclude oligos with too many degenerate bases
  oligos <- oligos[oligos$degenerates <= max_degenerates, ]
  # Exclude oligos with too high degeneracy
  oligos <- oligos[oligos$degeneracy <= max_degeneracy, ]
  # Identify and exclude oligos that are duplicated
  unique_oligos <- match(oligos$majority, unique(oligos$majority))
  oligos <- oligos[unique_oligos, ]
  return(oligos)
}

#' Calculate GC content and tm of oligos
#'
#' @param oligos A tibble with oligos.
#'
#' @param gc_range The GC-content range of each oligo. Can range between
#' 0 and 1. The default is \code{c(0.45, 0.55)}.
#'
#' @param tm_range The Tm-range of each oligo. Can range between 20 and 90.
#' The default is \code{c(48, 70)}.
#'
#' @param conc_oligo The concentration of oligonucleotide in M,
#' ranging from 0.2e-07 M (20 nM) to 2e-06 M (2000 nM).
#' The default value is 5e-07 M (500 nM) (for Tm calculation)
#'
#' @param conc_na The sodium ion concentration in M, ranging
#' from 0.01 M to 1 M. The default value is 0.05 M (50 mM)
#' (for Tm calculation).
#'
#' @section Excluded oligos:
#' The function excludes oligos with
#' more than than three consecutive runs of the same dinucleotide
#' (e.g. 'tatatata'), oligos with more than four consecutive runs of the
#' same nucleotide, and oligos that are duplicated.
#'
#' @section Tm:
#' The melting temperature is calculated using the nearest-neigbour method.
#' The oligo concentration is set to 500 nM and the sodium ion concentration
#' is set to 50 mM.
#'
#' Assumptions for Tm calculation:
#'
#' Oligos are not expected to be self-complementary, so no symmetry
#' correction is done.
#'
#' We assume that the oligo concentration is much higher
#' than the target concentration.
#'
#' See references for table values and equations.
#'
#' #Warning:
#' GC-content and Tm are calculated based on the majority oligos, and
#' may thus be misleading for degenerate (iupac) oligos.
#'
#' @return A tibble (a data frame) with oligo candidates.
#'
#' @noRd
add_gc_tm <- function(
                      oligos,
                      gc_range = c(0.45, 0.55),
                      tm_range = c(48, 70),
                      conc_oligo = 5e-07,
                      conc_na = 0.05) {
  if (!(min(gc_range) >= 0 && max(gc_range) <= 1)) {
    stop("gc_range must be between 0 and 1, e.g. c(0.45, 0.65)", call. = FALSE)
  }
  if (!(min(tm_range) >= 20 && max(tm_range) <= 90)) {
    stop("tm_range must be between 20 and 90, e.g. c(55, 60)", call. = FALSE)
  }
  # Calculate GC content of all majority oligos
  gc_majority <- purrr::map_dbl(oligos$majority, ~ gc_content(.x))
  oligos <- tibble::add_column(oligos, gc_majority)
  # Exclude oligos with GC content outside the stated thresholds
  oligos <- oligos[oligos$gc_majority >= min(gc_range), ]
  oligos <- oligos[oligos$gc_majority <= max(gc_range), ]
  # Calculate Tm of all majority oligos
  tm_majority <- tm(oligos$majority, conc_oligo = conc_oligo, conc_na = conc_na)
  oligos <- tibble::add_column(oligos, tm_majority)
  # Exclude oligos with Tm outside the stated thresholds
  oligos <- oligos[oligos$tm_majority >= min(tm_range), ]
  oligos <- oligos[oligos$tm_majority <= max(tm_range), ]
  return(oligos)
}

#' Convert a DNA sequence to a regular expression
#'
#' \code{make_regex} converts a DNA sequence
#' to a regular expression for pattern matching.
#'
#' @param x A DNA sequence (a character vector of length one), e.g. 'cttgtr'.
#'
#' @details Valid bases for \code{x} are 'a', 'c', 'g', 't', 'r', 'y', 'm',
#' 'k', 's', 'w', n', 'h', 'd', 'v', 'b' and '-'
#'
#' @return A regular expression of x (e.g. '(c)(t)(t)(g)(t)(a|g)').
#'
#' @examples
#' make_regex("cttrng")
#' @noRd
make_regex <- function(x) {
  if (typeof(x) != "character" || length(x) != 1) {
    stop("x must be a character vector of length one", call. = FALSE)
  }
  if (grepl("[^acgtryswkmbdhvn-]", x)) {
    stop("x contains at least one invalid base.
      Valid bases are 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w',
      'n', 'h', 'd', 'v', 'b' and '-'",
      call. = FALSE
    )
  }
  x <- split_sequence(x)
  # Go through each base of the DNA sequence
  regx <- purrr::map(x, function(i) {
    # Check which bases the IUPAC base at position 'i' corresponds to
    all_bases <- unname(degenerate_lookup[i])
    all_bases <- unlist(strsplit(all_bases, split = ","))
    return(all_bases)
  })
  regx <- purrr::map(regx, ~ paste(.x, collapse = "|"))
  regx <- purrr::map(regx, ~ paste0("(", .x, ")"))
  regx <- unlist(paste(regx, collapse = ""))
  return(regx)
}

#' Check if oligos matches their targets
#'
#' \code{check_match} checks if oligos matches with their
#' intended target sequences.
#'
#' @param x A tibble with oligos
#'
#' @param y The intended target. An alignment of DNA sequences
#' (an object of class 'rprimer_alignment').
#'
#' @return A tibble (a data frame)
#' with columns describing the proportion of perfectly matching
#' sequences for each oligo, and a column named
#' 'match_report', which contains a matrix with information about
#' which sequences the oligo matches perfectly to.
#'
#' @noRd
check_match <- function(x, y) {
  if (!inherits(y, "rprimer_alignment")) {
    stop("An rprimer_alignment object is expected for y.", call. = FALSE)
  }
  majority <- ifelse(
    !is.na(x$majority), x$majority,
    purrr::map_chr(x$majority_rc, reverse_complement)
  )
  iupac <- ifelse(
    !is.na(x$iupac), x$iupac, purrr::map_chr(x$iupac_rc, reverse_complement)
  )
  iupac <- purrr::map_chr(iupac, make_regex)
  # Shorten the sequence names to only accession numbers
  names(y) <- purrr::map_chr(names(y), truncate_name)
  # Make a matrix that describes which sequences the oligo matches perfectly to
  match_matrix <- purrr::map(seq_len(nrow(x)), function(i) {
    match_majority <- grepl(majority[[i]], y)
    match_iupac <- grepl(iupac[[i]], y)
    match <- cbind(match_majority, match_iupac)
    colnames(match) <- c("match_majority", "match_iupac")
    rownames(match) <- names(y)
    return(match)
  })
  names(match_matrix) <- x$iupac
  # Calculate the match percentage for each oligo
  match_percentage <- purrr::map(match_matrix, colMeans)
  match_percentage <- do.call("rbind", match_percentage)
  colnames(match_percentage) <- c("pm_majority", "pm_iupac")
  match_percentage <- tibble::as_tibble(match_percentage)
  # Add this information to the rprimer_oligo object
  x <- dplyr::bind_cols(x, match_percentage)
  match_matrix <- tibble::tibble(match_matrix)
  x <- dplyr::bind_cols(x, match_matrix)
  return(x)
}

# Context: Calculate tm =======================================================

#' Split a DNA sequence into nearest neighbors
#'
#' \code{nn_split} splits an oligo sequence into nearest neighbors
#' (for calculation of deltaG, deltaH and Tm)
#'
#' @param x a DNA sequence with at least two bases, e.g. 'ctta'
#' (a character vector of length one).
#'
#' @return The nearest neighbors of x (a character vector).
#'
#' @example
#' nn_split('ccgtncg')
#'
#' @noRd
nn_split <- function(x) {
  if (!(!is.na(x) && is.character(x) && nchar(x) > 1)) {
    stop("x must be a character vector of length one,
      with at least two characters (e.g. 'caaggnt')", call. = FALSE)
  }
  x <- split_sequence(x)
  from <- (seq_along(x) - 1)[-1]
  to <- seq_along(x)[-1]
  nn <- purrr::map_chr(from, ~ paste(x[from[[.x]]:to[[.x]]], collapse = ""))
  return(nn)
}

#' Calculate dH or dS of nearest neighbors using lookup tables
#'
#' @param x A character vector with nearest-neighbor pairs
#' of DNA sequences (e.g. \code{c('ct', 'tt', 'ta')})
#'
#' @param table The lookup table that should be used.
#' Either 'dH' (entropy) or 'dS' (enthalpy).
#'
#' @details \code{x} cannot contain other characters than
#' 'a', 'c', 'g', 't' and '-'.
#'
#' @return The corresponding values for dH or dS (in cal/M)
#'
#' @examples
#' nn_lookup(c("ac", "cc"), table = "dS")
#' @noRd
nn_lookup <- function(x, table) {
  valid_tables <- c("dH", "dS")
  if (!is.character(table) || !table %in% valid_tables) {
    stop("table must be set to either 'dH' or 'dS'", call. = FALSE)
  }
  if (any(grepl("[^acgt]", x))) {
    stop(
      "x contain at least one invalid base.
          Valid bases are a, c, g and t.",
      call. = FALSE
    )
  }
  if (table == "dH") {
    selected_table <- nearest_neighbor_lookup$dH
  } else {
    selected_table <- nearest_neighbor_lookup$dS
  }
  matching <- selected_table[match(x, nearest_neighbor_lookup$bases)]
  if (is.null(ncol(x))) {
    result <- matching
  } else {
    result <- matrix(matching, ncol = ncol(x), byrow = FALSE)
  }
  return(result)
}

#' Initiation of DNA sequences for Tm calculation
#'
#' @param x One or more DNA sequences (a character vector).
#'
#' @details \code{x} cannot contain other characters than
#' 'a', 'c', 'g', 't' and '-'.
#'
#' @return The initiaion values for x.
#'
#' @noRd
init_3end <- function(x) {
  if (grepl("(t|a)$", x)) {
    c(H = 2.3 * 1000, S = 4.1)
  } else {
    c(H = 0.1 * 1000, S = -2.8)
  }
}
init_5end <- function(x) {
  if (grepl("^(t|a)", x)) {
    c(H = 2.3 * 1000, S = 4.1)
  } else {
    c(H = 0.1 * 1000, S = -2.8)
  }
}

#' Melting temperature
#'
#' \code{tm} calculates the melting temperature of one or
#' more perfectly matching DNA duplexes (i.e. oligo-target duplexes),
#' using the nearest neigbour method.
#'
#' @param x One or more DNA sequences (a character vector).
#'
#' @param conc_oligo The concentration of oligonucleotide in M,
#' ranging from 0.2e-07 M (20 nM) to 2.0e-06 M (2000 nM).
#' The default value is 5e-07 M (500 nM).
#'
#' @param conc_na The sodium ion concentration in M, ranging
#' from 0.01 M to 1 M. The default value is 0.05 M (50 mM).
#'
#' @details \code{x} cannot contain other characters than
#' 'a', 'c', 'g', 't' and '-'. All oligos must be of equal length.
#'
#' @section Calculation of Tm:
#'
#' No symmetry correction is done
#' (oligos are not expected to be self-complementary).
#'
#' We assume that the oligo concentration is much higher
#' than the target concentration.
#'
#' @return The melting temperature(s) of x.
#'
#' @references
#'
#' SantaLucia, J, et al. (1996)
#' Improved Nearest-Neighbor Parameters for Predicting DNA Duplex Stability.
#' Biochemistry, 35: 3555-3562 (Formula and salt correction are from here)
#'
#' Allawi, H. & SantaLucia, J. (1997)
#' Thermodynamics and NMR of Internal G·T Mismatches in DNA.
#' Biochemistry, 34: 10581?\200?10594 ##############################################################
#' (Duplex initiation parameters are from here)
#'
#' SantaLucia, J (1998) A unified view of polymer,
#' dumbell, and oligonucleotide DNA nearest-neighbor thermodynamics.
#' Proc. Natl. Acad. Sci. USA, 95: 1460-1465. (Table values are from here)
#'
#' @examples
#' tm("acggtgcctac")
#' tm(c("acggtgcctac", "acggtggctgc"))
#' @noRd
tm <- function(x, conc_oligo = 5e-07, conc_na = 0.05) {
  if (!is.double(conc_oligo) || conc_oligo < 2e-07 || conc_oligo > 2.0e-06) {
    stop("The oligo concentration must be between
           0.2e-07 M (20 nM) and 2e-06 M (2000 nM)", call. = FALSE)
  }
  if (!is.double(conc_na) || conc_na < 0.01 || conc_na > 1) {
    stop("The Na+ concentration must be between 0.01 and 1 M", call. = FALSE)
  }
  x <- tolower(x)
  # Find initiation values
  init_H <- purrr::map_dbl(x, ~ init_5end(.x)[["H"]] + init_3end(.x)[["H"]])
  init_S <- purrr::map_dbl(x, ~ init_5end(.x)[["S"]] + init_3end(.x)[["S"]])
  nn <- purrr::map(x, nn_split)
  # Check oligo length
  oligo_length <- purrr::map_int(nn, length)
  # I made a matrix based tm-calculation,
  # which means that all oligos must be of the same length
  if (length(unique(oligo_length)) != 1) {
    stop("All oligos must be of equal length", call. = FALSE)
  }
  # Find nearest neighbor values for dH and dS
  nn <- do.call("rbind", nn)
  dH_result <- nn_lookup(nn, "dH")
  dS_result <- nn_lookup(nn, "dS")
  # Sum dH and dS
  sumdH <- rowSums(dH_result) + init_H
  sumdS <- rowSums(dS_result) + init_S
  # Correct delta S for salt conc.
  N <- nchar(x[[1]]) - 1 # Number of phosphates
  sumdS <- sumdS + 0.368 * N * log(conc_na)
  tm <- sumdH / (sumdS + gas_constant * log(conc_oligo))
  tm <- tm - 273.15
  return(tm)
}

# Context: Get assays =========================================================

#' Combine match matrices
#'
#' @param x A tibble with assays
#'
#' @return A tibble with assays, with information on perfect matches    #### Not tested!!!!
#' for fwd and rev primers combined.
#'
#' @noRd
combine_match_matrices <- function(x) {
  fwd <- x$match_matrix_fwd
  rev <- x$match_matrix_rev
  fwd <- purrr::map(fwd, function(x) {
    colnames(x) <- paste0(colnames(x), "_fwd")
    x
  })
  rev <- purrr::map(rev, function(x) {
    colnames(x) <- paste0(colnames(x), "_rev")
    x
  })
  match_matrix <- purrr::map(seq_len(nrow(x)), function(y) {
    match <- cbind(fwd[[y]], rev[[y]])
    majority_all <- rowSums(match[, c(1, 3)])
    iupac_all <- rowSums(match[, c(2, 4)])
    majority_all <- ifelse(majority_all == 2, TRUE, FALSE)
    iupac_all <- ifelse(iupac_all == 2, TRUE, FALSE)
    match <- cbind(match, majority_all, iupac_all)
    return(match)
  })
  match_percentage <- purrr::map(match_matrix, ~ colMeans(.x[, 5:6]))
  match_percentage <- do.call("rbind", match_percentage)
  match_percentage <- tibble::as_tibble(match_percentage)
  names(match_percentage) <- paste0("pm_", names(match_percentage))
  x <- dplyr::bind_cols(x, match_percentage)
  match_matrix <- tibble::tibble(match_matrix)
  x <- dplyr::bind_cols(x, match_matrix)
  drop <- c("match_matrix_fwd", "match_matrix_rev")
  x <- x[, !(names(x) %in% drop)]
  return(x)
}

# Context: Add probes =========================================================

#' Add probes to match matrices
#'
#' @param x A tibble with assays (with probes)
#'
#' @return A tibble with assays, with information on perfect matches    #### Not tested!!!! ######### PROBLEM!
#' for fwd and rev primers and probes combined.
#'
#' @noRd
add_probe_to_match_matrix <- function(x) {
  assays <- x$match_matrix
  probe <- x$match_matrix_pr
  probe <- purrr::map(probe, function(x) {
    colnames(x) <- paste0(colnames(x), "_pr")
    x
  })
  match_matrix <- purrr::map(seq_len(nrow(x)), function(y) {
    match <- cbind(assays[[y]], probe[[y]])
    match <- match[, -(5:6)]
    majority_all <- rowSums(match[, c(1, 3, 5)])
    iupac_all <- rowSums(match[, c(2, 4, 6)])
    majority_all <- ifelse(majority_all == 3, TRUE, FALSE)
    iupac_all <- ifelse(iupac_all == 3, TRUE, FALSE)
    match <- cbind(match, majority_all, iupac_all)
    return(match)
  })
  match_percentage <- purrr::map(match_matrix, ~ colMeans(.x[, 7:8]))
  match_percentage <- do.call("rbind", match_percentage)
  match_percentage <- tibble::as_tibble(match_percentage)
  names(match_percentage) <- paste0("pm_", names(match_percentage))
  drop <- c("match_matrix_pr", "match_matrix", "pm_majority_all", "pm_iupac_all")
  x <- x[, !(names(x) %in% drop)]
  x <- dplyr::bind_cols(x, match_percentage)
  match_matrix <- tibble::tibble(match_matrix)
  x <- dplyr::bind_cols(x, match_matrix)
  return(x)
}

# Context: Plots ==============================================================

#' Calculate running average
#'
#' \code{running_average} calculates the centered running average of a
#' numeric vector.
#'
#' @param x A numeric vector.
#'
#' @param size The number of observations in each average.
#' If \code{NULL}, the size will be set to the nearest positive, nonzero
#' integer to \code{length(x)/100}.
#'
#' @return A tibble with position and running average of x.
#'
#' @examples
#' running_average(c(1, 3, 2, 1, 3, 4, 5), size = 2)
#' running_average(c(10, 22.4, 14.2, 44, 32))
#' @noRd
running_average <- function(x, size = NULL) {
  if (is.null(size)) {
    size <- round(length(x) / 100)
    if (size == 0) {
      size <- 1
    }
  }
  stopifnot(is.numeric(x) && size <= length(x) && size >= 1)
  sums <- c(0, cumsum(x))
  average <- (sums[((size + 1):length(sums))] - sums[(1:(length(sums) - size))]) / size
  # First, we set the position to get a trailed moving average
  position <- (1 + size - 1):length(x)
  # Then, to get a centered running average,
  # we align each moving average to the midpoint of the range of observations
  # that each average includes
  midpoint <- size / 2
  position <- position - midpoint
  df <- tibble::tibble(position, average)
  return(df)
}

#' Calculate running average of GC content
#'
#' \code{gc_running_average} calculates the centered running average of
#' the GC-content of a DNA sequence
#'
#' @param x a DNA sequence (a character vector of length one).
#'
#' @param size The number of observations in each average.
#' If \code{NULL}, the size will be set to the nearest positive, nonzero
#' integer to \code{length(x)/100}.
#'
#' @details \code{x} cannot contain other characters than
#' 'a', 'c', 'g', 't' and '-'.
#'
#' @return A tibble with position and running average of x.
#'
#' @examples
#' running_average(c("gtttgtttctt"), size = 2)
#' running_average(c("cttggggtttctttggtt-ttagg"))
#' @seealso
#' gc_content
#'
#' @noRd
gc_running_average <- function(x, size = NULL) {
  if (is.null(size)) {
    size <- round(length(x) / 100)
    if (size == 0) {
      size <- 1
    }
  }
  stopifnot(is.character(x) && size <= length(x) && size >= 1)
  begins <- 1:(length(x) - size + 1)
  ends <- begins + size - 1
  frame <- 1:(length(x) - size + 1)
  average <- purrr::map_dbl(frame, function(y) {
    string <- paste(x[begins[[y]]:ends[[y]]], collapse = "")
    return(gc_content(string))
  })
  # First, we set the position to get a trailed moving average
  position <- (1 + size - 1):length(x)
  # Then, to get a centered running average,
  # we align each moving average to the midpoint of the range of observations
  # that each average includes
  midpoint <- size / 2
  position <- position - midpoint
  df <- tibble::tibble(position, average)
  return(df)
}

#' Draw a rectangle
#'
#' @param from Where (at the x-axis) the rectangle starts.
#'
#' @param to Where (at the x-axis) the rectangle ends.
#'
#' @return A rectangle.
#'
#' @noRd
rectangle <- function(from, to) {
  graphics::rect(
    from, 0, to, 5,
    border = NA,
    col = grDevices::rgb(123, 149, 169, alpha = 100, maxColorValue = 200), xpd = NA
  )
}

#' Sequence detail plot
#'
#' @param x An object of class 'rprimer_properties'
#'
#' @return A visual representation of \code{x}
#'
#' @noRd
sequence_detail_plot <- function(x) {
  identity_plot <- graphics::plot(
    x$position, x$identity,
    type = "h", ylim = c(0, 1),
    ylab = "identity", xlab = "", xaxt = "n",
    col = ifelse(x$identity < 1, "gray80", "gray60")
  )
  identity_line <- graphics::lines(running_average(x$identity))
  entropy_plot <- graphics::plot(
    x$position, x$entropy,
    type = "h",
    ylim = c(0, max(x$entropy, na.rm = TRUE) * 1.1),
    ylab = "shannon entropy",
    xlab = "", xaxt = "n", col = ifelse(x$entropy > 0, "gray80", "gray60")
  )
  entropy_line <- graphics::lines(running_average(x$entropy))
  gc_plot <- graphics::plot(
    x$position,
    pch = NA, ylab = "gc content", ylim = c(0, 1),
    xlab = "", xaxt = "n"
  )
  graphics::clip(0, nrow(x), -1, 2)
  graphics::abline(h = 0.5, col = "gray80")
  gc_line <- graphics::lines(gc_running_average(x$majority))
  gap_plot <- graphics::plot(
    x$position, x$gaps,
    type = "h", ylim = c(0, 1), ylab = "gaps", xlab = ""
  )
  identity_plot
  identity_line
  entropy_plot
  entropy_line
  gap_plot
  graphics::mtext(
    side = 1, outer = TRUE, "position in consensus sequence",
    line = 3, cex = 0.7
  )
}

#' Sequence barplot
#'
#' @param x A selection of a sequence profile
#'
#' @param ... Additional arguments that should be
#' passed to the barplot, e.g. a title.
#'
#' @return A barplot.
#'
#' @noRd
sequence_barplot <- function(x, ...) {
  colors1 <- c("#7B95A9", "#E7D0D8", "#D09F99", "#404038")
  names(colors1) <- c("a", "c", "g", "t")
  colors2 <- grDevices::gray.colors(nrow(x) - 4, start = 0.6)
  names(colors2) <- setdiff(rownames(x), names(colors1))
  colors <- c(colors1, colors2)
  # Make the plot
  graphics::barplot(
    x,
    space = 0, xaxt = "n", font.main = 1, border = "grey80",
    col = colors[rownames(x)], legend = TRUE, ylab = "Proportion",
    ylim = c(0, 1), ...,
    args.legend = list(
      x = "right", box.col = NA, border = "grey80",
      bg = grDevices::rgb(60, 60, 60,
        alpha = 50,
        maxColorValue = 200
      ), inset = c(0, -0.3)
    )
  )
}
