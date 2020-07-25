# Functions for the package rprimer (not exported) The code was written by Sofia in R 4.0.0

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

#' 'rprimer_alignment' constructor and validator
#'
#' \code{new_rprimer_alignment} creates an object of class
#' 'rprimer_alignment'
#'
#' @param x An rprimer_alignment-like object.
#'
#' @details an rprimer_alignment object must must be a list and
#' contain at least one DNA
#' sequence, and all sequences must be a character vector of lentgh one.
#' All sequences (including gaps) must be of the same length,
#' and the alignment must beat least 200 bases long.
#' All sequences must have unique names. Valid bases are
#' a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w',
#' 'n', 'h', 'd', 'v', 'b' and '-')
#'
#' @return An rprimer_alignment object if the validation is succeeds.
#' An error message if not.
#'
#' @noRd
new_rprimer_alignment <- function(x = list()) {
    # x must be a list
    stopifnot(is.list(x))
    # All sequences must be a character vector of length one
    has_length_one <- purrr::map_lgl(x, ~is.character(.x) && length(.x) == 1)
    if (any(has_length_one == FALSE)) {
        stop(
          "All sequences must be a character vector of length one",
          call. = FALSE
        )
    }
    # All sequences must contain the same number of characters
    sequence_lengths <- purrr::map_int(x, nchar)
    if (length(unique(sequence_lengths)) != 1) {
        stop(
        "The sequences does not appear to be aligned.
        All sequences (including gaps) are not of the same length.",
        call. = FALSE
      )
    }
    # The alignment must be at least 200 bases long
    if (unique(sequence_lengths) < 200) {
        stop(paste(
          "The alignment has", unique(sequence_lengths), "bases.
          The minumum is 200."
       ),
       call. = FALSE)
    }
    # All sequence names must be unique
    unique_name <- length(unique(names(x))) == length(x)
    if (unique_name == FALSE) {
        stop("All sequences must have unique names", call. = FALSE)
    }
    # All sequences must be in lowercase format and contain only valid bases
    non_valid_base <- purrr::map_lgl(x, ~grepl("[^acgtrymkswnhdvb-]", .x))
    if (any(non_valid_base)) {
        stop(
        "At least one sequence contain one or more invalid characters.
         Valid characters are 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w',
         'n', 'h', 'd', 'v', 'b' and '-')", call. = FALSE
      )
    }
    # Set class attribute
    x <- structure(x, class = "rprimer_alignment")
    return(x)
}

#' 'rprimer_sequence_profile' constructor and validator
#'
#' \code{new_rprimer_sequence_profile} creates an object of class
#' 'rprimer_sequence_profile'
#'
#' @param x An rprimer_sequence_profile-like object.
#'
#' @details an rprimer_sequence_profile object must be a matrix of type
#' 'double'. Each element must have a value between 0 and 1. The matrix
#' must have both rownames and colnames.
#'
#' @return An rprimer_sequence_profile object if the validation is succeeds.
#' An error message if not.
#'
#' @noRd
new_rprimer_sequence_profile <- function(x = matrix()) {
    # Integrity checks
    stopifnot(
      is.matrix(x), is.double(x), max(x) <= 1, min(x) >= 0,
      !is.null(rownames(x)), !is.null(colnames(x))
    )
    # Set class attribute
    x <- structure(x, class = "rprimer_sequence_profile")
    return(x)
}

#' Majority consensus sequence
#'
#' \code{majority_consensus} returns the majority consensus sequence of an
#' alignment of DNA sequences.
#'
#' @param x A nucleotide profile of class 'rprimer_sequence_profile', i.e.
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
    if (!inherits(x, "rprimer_sequence_profile")) {
        stop("An rprimer_sequence_profile object is expected.", call. = FALSE)
    }
    # Function to identify the most common base at a position
    find_most_common_base <- function(x, y) {
        most_common <- rownames(x)[which(y == max(y))]
        # If there are ties, the most common base will be randomly selected
        if (length(most_common > 1))
            most_common <- sample(most_common, 1)
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
#' as_iupac('a,c')
#' as_iupac('r')
#' as_iupac('tg') # Will return NA since the bases are not separated by comma
#'
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
    bases <- c("a", "c", "g", "t", "-")  # accepted nucleotides
    match <- x %in% bases
    x <- x[!(match == FALSE)]  # exclude non accepted nucleotides
    x <- paste(x, collapse = ",")
    iupac <- unname(iupac_lookup[x]) # match with lookup table
    return(iupac)
}

#' Iupac consensus sequence
#'
#' \code{iupac_consensus} returns the iupac consensus sequence from an
#' alignment of DNA sequences.
#'
#' @param x A nucleotide profile of class 'rprimer_sequence_profile', i.e.
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
    if (!inherits(x, "rprimer_sequence_profile")) {
        stop("An rprimer_sequence_profile object is expected.", call. = FALSE)
    }
    if (!is.double(threshold) || threshold < 0 || threshold > 0.2) {
        stop(paste0(
        "The threshold was set to ", threshold, ". Valid threshold values
     are from 0 to 0.2"
        ))
    }
    # Select only a, c, g, t and - to count in the iupac consensus sequence
    bases <- c("a", "c", "g", "t", "-")
    x <- x[which(rownames(x) %in% bases), ]
    bases_to_include <- apply(x, 2, function(y) {
      paste(rownames(x)[which(y >= threshold)], collapse = ",")
    })
    bases_to_include <- unname(bases_to_include)
    consensus <- purrr::map_chr(bases_to_include, ~as_iupac(.x))
    if (any(is.na(consensus))) {
        warning("The consensus sequence contain NAs.
    Try to lower the threshold value.", call. = FALSE)  # Check if this works!
    }
    return(consensus)
}

#' Gap frequency
#'
#' \code{gap_frequency} returns the gap frequency at each position in an
#' alignment of DNA sequences.
#'
#' @param x A nucleotide profile of class 'rprimer_sequence_profile', i.e.
#' a numeric m x n matrix, where m is the number of nucleotides in the
#' alignment, and n is the number of positions in the alignment.
#'
#' @details Gaps are recognised as "-".
#'
#' @return the gap frequency (a numeric vector of length n).
#'
#' @noRd
gap_frequency <- function(x) {
    if (!inherits(x, "rprimer_sequence_profile")) {
        stop("An rprimer_sequence_profile object is expected.", call. = FALSE)
    }
    gaps <- x[which(rownames(x) == "-"), ]
    gaps <- unname(gaps)
    return(gaps)
}

#' Nucleotide identity
#'
#' #' \code{nucleotide_identity} returns the nucleotide identity from an
#' alignment of DNA sequences.
#'
#' @param x A nucleotide profile of class 'rprimer_sequence_profile', i.e.
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
    if (!inherits(x, "rprimer_sequence_profile")) {
        stop("An rprimer_sequence_profile object is expected.", call. = FALSE)
    }
    # We want to assess identity based on DNA bases,
    # i.e. ignore gaps and degenerate positions, so we make a subset (s) of x
    # with the rows named a, c, g and t.
    bases <- c("a", "c", "g", "t")
    s <- x[which(rownames(x) %in% bases), ]
    # Calculate relative proportions of the bases in s
    s <- apply(s, 2, function(x) x / sum(x))
    # Find the largest proportion at each position
    identity <- apply(s, 2, max)
    identity <- unname(identity)
    identity[which(is.na(identity))] <- 0
    return(identity)
}

#' Shannon entropy
#'
#' #' \code{shannon_entropy} returns the shannon entropy from an
#' alignment of DNA sequences.
#'
#' @param x A nucleotide profile of class 'rprimer_sequence_profile', i.e.
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
    if (!inherits(x, "rprimer_sequence_profile")) {
        stop("An rprimer_sequence_profile object is expected.", call. = FALSE)
    }
    # We want to assess identity based on DNA bases,
    # i.e. ignore gaps and degenerate positions, so we make a subset (s) of x
    # with the rows named a, c, g and t.
    bases <- c("a", "c", "g", "t")
    s <- x[which(rownames(x) %in% bases), ]
    # Calculate relative proportions of the bases in s
    s <- apply(s, 2, function(x) x / sum(x))
    entropy <- apply(s[which(rownames(s) %in% bases), ], 2, function(x) {
        ifelse(x == 0, 0, x * log2(x))
    })
    entropy <- -colSums(entropy)
    entropy <- unname(entropy)
    entropy <- abs(entropy)  # abs to avoid -0 (due to neg sums)
    entropy[which(is.na(entropy))] <- 0
    return(entropy)
}

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
#' get_nmers(c('c', 'g', 't', 't', 'c', 'g'), n = 2)
#'
#' @noRd
get_nmers <- function(x, n = NULL) {
    if (!is.character(x)) {
        stop("x must be a character vector.", call. = FALSE)
    }
    if (is.null(n)) {
        n <- round(length(x) / 10)
        if (n == 0)
            n <- 1
    }
    if (!is.numeric(n) || n < 1 || n > length(x)) {
        stop("n must be an integer between 1 and length(x)", call. = FALSE)
    }
    begin <- 1:(length(x) - n + 1)
    end <- begin + n - 1
    nmer <- purrr::map_chr(
      begin, ~paste(x[begin[[.x]]:end[[.x]]], collapse = "")
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
#' gc_content('acgttcc')
#' gc_content('acgttcc--')
#' gc_content('acgrn') # Will return an error because of an invalid base.
#'
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
    x <- strsplit(x, split = "")
    x <- unlist(x, use.names = FALSE)
    gc_count <- length(which(x == "c" | x == "g"))
    # Gaps will not be included in the total count
    total_count <- length(which(x == "a" | x == "c" | x == "g" | x == "t"))
    gc <- gc_count / total_count
    return(gc)
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
#' @return The complement sequence of x. Non valid bases will return as NA.
#'
#' @examples
#' reverse_complement('cttgtr')
#'
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
            call. = FALSE)
    }
    x <- strsplit(x, split = "")
    complement <- complement_lookup[unlist(x)]
    complement <- unname(complement)
    return(complement)
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
#' reverse_complement('cttgtr')
#'
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
            call. = FALSE)
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
#'
#' @noRd
running_sum <- function(x, n = NULL) {
    if (!(is.numeric(x) && is.vector(x))) {
        stop("A numeric vector is expected for x.", call. = FALSE)
    }
    if (is.null(n)) {
        n <- round(length(x) / 10)
        if (n == 0)
            n <- 1
    }
    if (!is.numeric(n) || n < 1 || n > length(x)) {
        stop("n must be a number between 1 and length(x)", call. = FALSE)
    }
    cumul <- c(0, cumsum(x))
    runsum <- cumul[(n + 1):length(cumul)] - cumul[1:(length(cumul) - n)]
    return(runsum)
}

#' Regular expression to repeat a pattern
#'
#' \code{repeat_pattern} creates a regular expression for repeating a pattern.
#'
#' @param n An integer equal to or greater than 1.
#' The number of times the pattern is to be printed.
#'
#' @return A regular expression that can be pasted with another regular
#' expression, so that it can be repeated \code{n} times.
#'
#' @examples
#' repeat_pattern(3) # this pattern will be printed 3 times when pasted
#' to a regular expression.
#'
#' @noRd
repeat_pattern <- function(n) {
    if (!is.numeric(n) || n < 1) {
        stop("n must be a number, and at least 1.", call. = FALSE)
    }
    regex <- paste(rep("\\1", n - 1), collapse = "")
    return(regex)
}

#' Exclude oligos with a certain pattern
#'
#' \code{exclude_oligos} replaces oligos with a certain pattern
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
#' exclude_oligos(c('cttgttatttt', 'cgattctg'), "a$|t$")
#'
#' @noRd
exclude_oligos <- function(x, pattern) {
  if (!is.character(x)) {
    stop("x must be a character vector", call. = FALSE)
  }
  if (!is.character(pattern)) {
    stop("pattern must be a character vector", call. = FALSE)
  }
  regex <- pattern
  x[which(grepl(regex, x))] <- NA
  return(x)
}

#' Exclude unwanted oligos
#'
#' \code{exclude unwanted oligos} replaces oligos with many consecutive
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
#' exclude_unwanted_oligos(c('gtgtaaaaatatttt', 'cgattctg'))
#'
#' @noRd
exclude_unwanted_oligos <- function(
  x, avoid_3end_ta = FALSE, avoid_5end_g = FALSE, avoid_3end_runs = FALSE
  ) {
    # Remove oligos with at least 4 'runs' of the same dinucleotide
    x <- exclude_oligos(
      x, paste0("(at|ta|ac|ca|ag|ga|gt|tg|cg|gc|tc|ct)", repeat_pattern(4))
    )
    # Remove oligos with at least 5 'runs' of the same nucleotide
    x <- exclude_oligos(x, paste0("([a-z])", repeat_pattern(5)))
    if (avoid_3end_ta == TRUE) {
        # Remove oligos with t or a at the 3' end
        x <- exclude_oligos(x,"a$|t$")
    }
    if (avoid_5end_g == TRUE) {
        # Remove oligos with g at the 5' end
        x <- exclude_oligos(x, "^g")
    }
    if (avoid_3end_runs == TRUE) {
        # Remove oligos with at least 3 'runs' of the same nucleotide in 3'
        x <- exclude_oligos(x, paste0("([a-z])", repeat_pattern(3), "$"))
    }
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
#' @return the number of degenerate bases in x (an integer).
#'
#' @examples
#' count_degenerates('cttnra')
#'
#' @noRd
count_degenerates <- function(x) {
    if (typeof(x) != "character" || length(x) != 1) {
        stop("x must be a character vector of length one", call. = FALSE)
    }
    if (grepl("[^acgtrymkswnhdvb-]", x)) {
        stop("x contains at least one invalid base.
      Valid bases are 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w',
      'n', 'h', 'd', 'v', 'b' and '-'",
            call. = FALSE)
    }
    nt <- c("a", "c", "g", "t", "-")
    x <- strsplit(x, split = "")
    x <- unlist(x)
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
#' count_degeneracy('cttnra')
#'
#' @noRd
count_degeneracy <- function(x) {
    if (typeof(x) != "character" || length(x) != 1) {
        stop("x must be a character vector of length one", call. = FALSE)
    }
    if (grepl("[^acgtrymkswnhdvb-]", x)) {
        stop("x contains at least one invalid base.
      Valid bases are 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w',
      'n', 'h', 'd', 'v', 'b' and '-'",
            call. = FALSE)
    }
    x <- strsplit(x, split = "")
    x <- unlist(x, use.names = FALSE)
    # Find the number of nucleotides at each position in x
    n_nucleotides <- degeneracy_lookup[x]
    # Calculate the total number of DNA sequences in x
    degeneracy <- prod(n_nucleotides)
    return(degeneracy)
}

#' Split a DNA sequence into nearest neighbours
#'
#' \code{nn_split} splits an oligo sequence into nearest neighbours
#' (for calculation of deltaG, deltaH and Tm)
#'
#' @param x a DNA sequence with at least two bases, e.g. 'ctta'
#' (a character vector of length one).
#'
#' @return The nearest neighbours of x (a character vector).
#'
#' @example
#' nn_split('ccgtncg')
#'
#' @noRd
nn_split <- function(x) {
    if (!(!is.na(x) && is.character(x) && is.vector(x) && nchar(x) > 1)) {
        stop("x must be a character vector of length one,
      with at least two characters (e.g. 'caaggnt')", call. = FALSE)
    }
    x <- strsplit(x, split = "")
    x <- unlist(x, use.names = FALSE)
    from <- 1:(length(x) - 1)
    to <- 2:length(x)
    nn <- purrr::map_chr(from, ~paste(x[from[[.x]]:to[[.x]]], collapse = ""))
    return(nn)
}

#' Calculate dH or dS of nearest neighbours using lookup tables
#'
#' @param x A character vector with nearest-neigbour pairs
#' of DNA sequences (e.g. \code{c('ct', 'tt', 'ta')})
#'
#' @param table The lookup table that should be used.
#' Either 'dH' (entropy) or 'dS' (enthalpy).
#'
#' @details \code{x} cannot contain other characters than
#' 'a', 'c', 'g', 't' and '-'.
#'
#' @return The corresponding values for dH or dS,
#'
#' @examples
#' nn_lookup(c('ac', 'cc'), table = 'dS')
#'
#' @noRd
nn_lookup <- function(x, table) {
    valid_tables <- c("dH", "dS")
    if (!is.character(table) || !table %in% valid_tables) {
        stop("table must be set to either 'dH' or 'dS'", call. = FALSE)
    }
    if (any(grepl("[^acgt.]", x))) {
        stop("x contain at least one invalid base. Valid bases are a, c, g and t.", call. = FALSE)
    }
    if (table == "dH") {
        selected_table <- nearest_neighbour_lookup$dH
    } else {
        selected_table <- nearest_neighbour_lookup$dS
    }
    matching <- selected_table[match(x, nearest_neighbour_lookup$bases)]
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
    if (grepl("^(t|a)", x))
        c(H = 2.3 * 1000, S = 4.1) else c(H = 0.1 * 1000, S = -2.8)
}
init_5end <- function(x) {
    if (grepl("(t|a)$", x))
        c(H = 2.3 * 1000, S = 4.1) else c(H = 0.1 * 1000, S = -2.8)
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
#' ranging from 0.2e-07 M (20 nM) to 2e-06 M (2000 nM).
#' The default value is 5e-07 M (500 nM).
#'
#' @param conc_na The sodium ion concentration in M, ranging
#' from 0.01 M to 1 M. The default value is 0.05 M (50 mM).
#'
#' @details \code{x} cannot contain other characters than
#' 'a', 'c', 'g', 't' and '-'.
#' All oligos in x must be of equal length.
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
#' Owczary et al. (2004)
#' Effects of Sodium Ions on DNA Duplex Oligomers:?\200? # NOT IN USE !
#' Improved Predictions of Melting Temperatures.
#' Biochemistry 43: 3537-3554
#'
#' SantaLucia, J, et al. (1996)
#' Improved Nearest-Neighbor Parameters for Predicting DNA Duplex Stability.
#' Biochemistry, 35: 3555-3562 (Formula and salt correction are from here)
#'
#' Allawi, H. & SantaLucia, J. (1997)
#' Thermodynamics and NMR of Internal G·T Mismatches in DNA.
#' Biochemistry, 34: 10581?\200?10594 (Duplex initiation parameters are from here)
#'
#' #' SantaLucia, J (1998) A unified view of polymer,
#' dumbell, and oligonucleotide DNA nearest-neighbour thermodynamics.
#' Proc. Natl. Acad. Sci. USA, 95: 1460-1465. (Table values are from here)
#'
#' @examples
#' tm('acggtgcctac')
#' tm(c('acggtgcctac', 'acggtggctgc'))
#'
#' @noRd
tm <- function(x, conc_oligo = 5e-07, conc_na = 0.05) {
    init_H <- purrr::map_dbl(x, ~init_5end(.x)[["H"]] + init_3end(.x)[["H"]])
    init_S <- purrr::map_dbl(x, ~init_5end(.x)[["S"]] + init_3end(.x)[["S"]])
    nn <- purrr::map(x, nn_split)
    # Check oligo length
    oligo_length <- purrr::map_int(nn, length)
    # I made a matrix based tm-calculation, which means that all oligos must be of the same length
    if (length(unique(oligo_length)) != 1) {
        stop("All oligos are not of the same length", call. = FALSE)
    }
    nn <- do.call("rbind", nn)
    dH_result <- nn_lookup(nn, "dH")
    dS_result <- nn_lookup(nn, "dS")
    sumdH <- rowSums(dH_result) + init_H
    sumdS <- rowSums(dS_result) + init_S
    # Correct delta S for salt conc
    N <- nchar(x[[1]]) - 1  # Number of phosphates
    sumdS <- sumdS + 0.368 * N * log(conc_na)
    tm <- sumdH / (sumdS + gas_constant * log(conc_oligo))
    tm <- tm - 273.15
    return(tm)
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
#' make_regex('cttrng')
#'
#' @noRd
make_regex <- function(x) {
    if (typeof(x) != "character" || length(x) != 1) {
        stop("x must be a character vector of length one", call. = FALSE)
    }
    if (grepl("[^acgtryswkmbdhvn-]", x)) {
        stop("x contains at least one invalid base.
      Valid bases are 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w',
      'n', 'h', 'd', 'v', 'b' and '-'",
            call. = FALSE)
    }
    x <- split_sequence(x)
    # Go through each base of the DNA sequence
    regx <- purrr::map(x, function(i) {
        # Check which bases the IUPAC base at position i correspond to by using a lookup table
        all_bases <- unname(degenerate_lookup[i])
        all_bases <- unlist(strsplit(all_bases, split = ","))
        return(all_bases)
    })
    regx <- purrr::map(regx, ~paste(.x, collapse = "|"))
    regx <- purrr::map(regx, ~paste0("(", .x, ")"))
    regx <- unlist(paste(regx, collapse = ""))
    return(regx)
}

#' Get all combinations of a DNA sequence with degenerate bases
#'
#' \code{expand_degenerate} returns all combinations of a DNA sequence with
#' degenerate bases.
#'
#' @param x A DNA sequence (a character vector of length one), e.g. 'cttgg'.
#' Valid bases are 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w',
#' n', 'h', 'd', 'v', 'b' and '-'.
#'
#' @return A list containing all combinations of the DNA sequence with
#'degenerate bases.
#'
#' @examples
#' expand_degenerate('cgtrn')
#'
#' @noRd
expand_degenerate <- function(x) {
    if (typeof(x) != "character" || length(x) != 1) {
        stop("x must be a character vector of length one", call. = FALSE)
    }
    if (grepl("[^acgtrymkswnhdvb-]", x)) {
        stop("x contains at least one invalid base.
      Valid bases are 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w',
      'n', 'h', 'd', 'v', 'b' and '-'",
            call. = FALSE)
    }
    x <- split_sequence(x)
    # Go through each base of the DNA sequence
    expanded <- purrr::map(x, function(i) {
        # Check which bases the IUPAC base at position i correspond to by using a lookup table
        all_bases <- unname(degenerate_lookup[[i]])
        all_bases <- unlist(strsplit(all_bases, split = ","))
        return(all_bases)
    })
    # Get all possible combinations of DNA sequences
    expanded <- expand.grid(expanded[seq_along(expanded)], stringsAsFactors = FALSE)
    expanded <- purrr::map(seq_len(nrow(expanded)), ~paste(expanded[.x, ], collapse = ""))
    expanded <- unlist(expanded, use.names = FALSE)
    return(expanded)
}

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
#' running_average(c(1,3,2,1,3,4,5), size = 2)
#' running_average(c(10, 22.4, 14.2, 44, 32))
#'
#' @noRd
running_average <- function(x, size = NULL) {
    if (is.null(size)) {
        size <- round(length(x)/100)
        if (size == 0)
            size <- 1
    }
    stopifnot(is.numeric(x) && size <= length(x) && size >= 1 && length(size == 1))
    sums <- c(0, cumsum(x))
    average <- (sums[((size + 1):length(sums))] - sums[(1:(length(sums) - size))])/size
    # First, we set the position to get a trailed moving average
    position <- (1 + size - 1):length(x)
    # Then, to get a centered running average, we align each moving average to the midpoint of the range of observations
    # that each average includes
    midpoint <- size/2
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
#' running_average(c('gtttgtttctt'), size = 2)
#' running_average(c('cttggggtttctttggtt-ttagg'))
#'
#' @seealso
#' gc_content
#'
#' @noRd
gc_running_average <- function(x, size = NULL) {
    if (is.null(size)) {
        size <- round(length(x)/100)
        if (size == 0)
            size <- 1
    }
    stopifnot(is.character(x) && size <= length(x) && size >= 1 && length(size == 1))
    begins <- 1:(length(x) - size + 1)
    ends <- begins + size - 1
    frame <- 1:(length(x) - size + 1)
    average <- purrr::map_dbl(frame, function(y) {
        string <- paste(x[begins[[y]]:ends[[y]]], collapse = "")
        return(gc_content(string))
    })
    # First, we set the position to get a trailed moving average
    position <- (1 + size - 1):length(x)
    # Then, to get a centered running average, we align each moving average to the midpoint of the range of observations
    # that each average includes
    midpoint <- size/2
    position <- position - midpoint
    df <- tibble::tibble(position, average)
    return(df)
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

#' Check if primers within an assay match their targets
#'
#' \code{check_primer_match} checks if primers in assays matches with their
#' intended target sequences.
#'
#' @param x Probes within an object of class
#' 'rprimer_oligo' or 'rprimer_assay'.
#'
#' @param y The intended target. An alignment of DNA sequences
#' (an object of class 'rprimer_alignment').
#'
#' @return A match report for each primer pair
#' (a list of logical matrices).
#'
#' @noRd
check_primer_match <- function(x, y) {
    # Collect the primer sequences from x
    primers <- x[c("majority_fwd", "iupac_fwd", "majority_rev", "iupac_rev")]
    # Reverse complement the reverse primers so that they can be matched against their targets
    primers$majority_rev <- purrr::map_chr(primers$majority_rev, reverse_complement)
    primers$iupac_rev <- purrr::map_chr(primers$iupac_rev, reverse_complement)
    # Select the unique primer pairs (will speed up the process a bit)
    primers <- unique(primers)
    # Make a matrix and convert primer sequences to regular expressions
    primer_matrix <- as.matrix(primers)
    primer_matrix <- apply(primer_matrix, c(1, 2), make_regex)
    # For each primer, check if it matches its targets
    match_matrix <- apply(primer_matrix, 1, function(x) {
        match <- purrr::map(x, ~grepl(.x, y))
        return(match)
    })
    match_matrix <- purrr::map(match_matrix, ~do.call("cbind", .x))
    # Set the sequence names as row names
    match_matrix <- purrr::map(match_matrix, function(x) {
        rownames(x) <- names(y)
        x
    })
    # Name the list with the primer pair so that it can be matched against x
    primers$majority_rev <- purrr::map_chr(primers$majority_rev, reverse_complement)
    names(match_matrix) <- paste0(primers$majority_fwd, "_", primers$majority_rev)
    # Expand to match x
    to_expand <- match(paste0(x$majority_fwd, "_", x$majority_rev), names(match_matrix))
    match_matrix <- match_matrix[to_expand]
    return(match_matrix)
}

#' Check if probes within an assay match their targets
#'
#' \code{check_probe_match} checks if probes in assays matches with their
#' intended target sequences.
#'
#' @param x Probes within an object of class
#' 'rprimer_oligo' or 'rprimer_assay'.
#'
#' @param y The intended target. An alignment of DNA sequences
#' (an object of class 'rprimer_alignment').
#'
#' @return A match report for each probe (a list of logical matrices).
#'
#' @noRd
check_probe_match <- function(x, y) {
    probes <- x[c("majority_pr", "iupac_pr", "sense_pr")]
    probes$majority_pr[which(probes$sense_pr == "neg")] <- purrr::map_chr(probes$majority_pr[which(probes$sense_pr == "neg")],
        reverse_complement)
    probes$iupac_pr[which(probes$sense_pr == "neg")] <- purrr::map_chr(probes$iupac_pr[which(probes$sense_pr == "neg")],
        reverse_complement)
    probes <- as.matrix(probes[, -ncol(probes)])
    probes <- apply(probes, c(1, 2), make_regex)
    match_matrix <- apply(probes, 1, function(x) {
        match <- purrr::map(x, ~grepl(.x, y))
        return(match)
    })
    match_matrix <- purrr::map(match_matrix, ~do.call("cbind", .x))
    match_matrix <- purrr::map(match_matrix, function(x) {
        rownames(x) <- names(y)
        x
    })
    return(match_matrix)
}

#' Draw a rectangle
#'
#'
#' @noRd
rectangle <- function(from, to) {
    rect(from, 0, to, 5, border = NA, col = rgb(123, 149, 169, alpha = 100, maxColorValue = 200), xpd = NA)
}

#' Make a sequence barplot
#'
#'
#' @param x A selection of a sequence profile
#'
#' @noRd
sequence_barplot <- function(x, ...) {
    colors1 <- c("#7B95A9", "#E7D0D8", "#D09F99", "#404038")
    names(colors1) <- c("a", "c", "g", "t")
    colors2 <- gray.colors(nrow(x) - 4, start = 0.6)
    names(colors2) <- setdiff(rownames(x), names(colors1))
    colors <- c(colors1, colors2)
    # Make the plot
    barplot(x, space = 0, xaxt = "n", font.main = 1, border = "grey80", col = colors[rownames(x)], legend = TRUE, ylab = "Proportion",
        ylim = c(0, 1), ..., args.legend = list(x = "right", box.col = NA, border = "grey80", bg = rgb(60, 60, 60, alpha = 50,
            maxColorValue = 200), inset = c(0, -0.3)))
}

#'  Round all doubles in a data frame
#'
#'  @param x A data frame/tibble.
#'
#'  @return A data frame/tibble where all vectors
#'  of type 'double' have been rounded to two digits.
#'
#'  @noRd
round_df_dbl <- function(x) {
    dbls <- purrr::map_lgl(x, is.double)
    dbls <- unname(dbls)
    dbls <- which(dbls == TRUE)
    x[dbls] <- round(x[dbls], 2)
    return(x)
}

#############################################################################

# Functions to create report
#' @noRd
print_assay_report <- function(x) {
    # x rprimer assay obj - one row
    if (any(grepl("match_matrix", names(x)))) {
        return(list(print_assay(x), print_degenerates(x), print_match_matrix(x$match_matrix[[1]])))
    } else {
        return(list(print_assay(x), print_degenerates(x)))
    }
}

#' @noRd
assay_detail_plot <- function(x, y) {
    # x = rprimer sequence profile y = assay object (one row)
    x <- x[which(rownames(x) != "-"), ]
    fwd <- x[, y$begin_fwd:y$end_fwd] #seq_len
    rev <- x[, y$begin_rev:y$end_rev] #seq-len
    rev <- rev[, ncol(rev):1]
    rownames(rev) <- unname(complement_lookup[rownames(rev)])
    if (any(grepl("_pr", names(y)))) {
        pr <- x[, y$begin_pr:y$end_pr]
        if (y$sense_pr == "neg") {
            pr <- pr[, ncol(pr):1]
            rownames(pr) <- unname(complement_lookup[rownames(pr)])
        }
    }
    if (any(grepl("_pr", names(y)))) {
        op <- par(mar = c(0.75, 4.57, 4.57, 0.75), mfrow = c(1, 3))
        on.exit(par(op))
        sequence_barplot(fwd, main = "forward")
        sequence_barplot(rev, main = "reverse")
        sequence_barplot(pr, main = "probe")
    } else {
        op <- par(mar = c(0.75, 4.57, 4.57, 0.75), mfrow = c(1, 2))
        on.exit(par(op))
        sequence_barplot(fwd, main = "forward")
        sequence_barplot(rev, main = "reverse")
    }
}

#' @noRd
assay_overview_plot <- function(x, y) {
    # x rprimer_sequence properties # y assay object (one row)
    op <- par(mfrow = c(4, 1), mai = c(0.1, 1, 0.1, 1), xpd = FALSE, oma = c(4, 0, 1, 0), mar = c(0.2, 4.1, 0.2, 2.1))
    on.exit(par(op))
    identity_plot <- plot(x$position, x$identity, type = "h", ylim = c(0, 1), ylab = "identity", xlab = "", xaxt = "n",
        col = ifelse(x$identity < 1, "gray80", "gray60"))
    identity_line <- lines(running_average(x$identity))
    entropy_plot <- plot(x$position, x$entropy, type = "h", ylim = c(0, max(x$entropy, na.rm = TRUE) * 1.1), ylab = "shannon entropy",
        xlab = "", xaxt = "n", col = ifelse(x$entropy > 0, "gray80", "gray60"))
    entropy_line <- lines(running_average(x$entropy))
    gc_plot <- plot(x$position, pch = NA, ylab = "gc content", ylim = c(0, 1), xlab = "", xaxt = "n")
    clip(0, nrow(x), -1, 2)
    abline(h = 0.5, col = "gray80")
    gc_line <- lines(gc_running_average(x$majority))
    gap_plot <- plot(x$position, x$gaps, type = "h", ylim = c(0, 1), ylab = "gaps", xlab = "")
    identity_plot
    identity_line
    entropy_plot
    entropy_line
    gap_plot
    rectangle(y$begin_fwd, y$end_fwd)
    rectangle(y$begin_rev, y$end_rev)
    if (any(grepl("_pr$", names(x))))
        rectangle(y$begin_pr, y$end_pr)
    mtext(side = 1, outer = TRUE, "position in consensus sequence", line = 3, cex = 0.7)
}

#' @noRd
print_assay <- function(x) {
    # x assay object (one row)
    if (any(grepl("match_matrix", names(x)))) {
        x <- x[which(!grepl("match_matrix", names(x)))]
    }
    all <- x[which(!grepl("_fwd$|_rev$|_pr$", names(x)))]
    names(all) <- gsub("_all$", "", names(all))
    fwd <- x[which(grepl("_fwd$", names(x)))]
    names(fwd) <- gsub("_fwd$", "", names(fwd))
    rev <- x[which(grepl("_rev$", names(x)))]
    names(rev) <- gsub("_rev$", "", names(rev))
    if (any(grepl("_pr", names(x)))) {
        pr <- x[which(grepl("_pr$", names(x)))]
        names(pr) <- gsub("_pr$", "", names(pr))
        my_list <- list(general = all, forward = fwd, reverse = rev, probe = pr)
    } else {
        my_list <- list(general = all, forward = fwd, reverse = rev)
    }
    return(my_list)
}

#' @noRd
print_degenerates <- function(x) {
    # x assay object (one row)
    iupac_oligos <- x[grep("^iupac", names(x))]
    names(iupac_oligos) <- gsub("iupac_", "", names(iupac_oligos))
    to_print <- purrr::map(iupac_oligos, function(x) {
        variant <- expand_degenerate(x)
        gc_content <- purrr::map_dbl(variant, ~gc_content(.x))
        tm <- tm(variant)
        table <- tibble::tibble(variant, gc_content, tm)
        mean_values <- round(c(gc_content = mean(gc_content), tm = mean(tm)), 2)
        return(list(all_variants = table, mean_values = mean_values))
    })
    return(to_print)
}

#' @noRd
print_match_matrix <- function(x) {
    # assay obj (one row)
    iupac <- x[, grep("iupac", colnames(x))]
    majority <- x[, grep("majority", colnames(x))]
    pm_majority <- majority[, which(colnames(majority) == "majority_all")]
    pm_majority <- names(pm_majority)[which(pm_majority == TRUE)]
    if (length(pm_majority) != 0) {
        to_remove_majority <- match(pm_majority, rownames(majority))
        mm_majority <- majority[-to_remove_majority, -ncol(majority)]
    } else {
        mm_majority <- majority
    }
    colnames(mm_majority) <- gsub("iupac", "", colnames(mm_majority))
    pm_iupac <- iupac[, which(colnames(iupac) == "iupac_all")]
    pm_iupac <- names(pm_iupac)[which(pm_iupac == TRUE)]
    if (length(pm_majority) != 0) {
        to_remove_iupac <- match(pm_iupac, rownames(iupac))
        mm_iupac <- iupac[-to_remove_iupac, -ncol(iupac)]
    } else {
        mm_iupac <- iupac
    }
    colnames(mm_iupac) <- gsub("majority", "", colnames(mm_iupac))
    result <- list(matches_perfectly_majority = pm_majority, mismatches_majority = mm_majority, matches_perfectly_iupac = pm_iupac,
        mismatches_iupac = mm_iupac)
    return(result)
}


