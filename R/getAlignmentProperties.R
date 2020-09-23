#' Get alignment properties
#'
#' \code{getAlignmentProperties} returns information from
#' an alignment profile
#'
#' @param x An alignment profile (an RprimerProfile object).
#'
#' @param IUPACThreshold
#' A number (0, 0.2].
#' At each position, all nucleotides with a proportion
#' \code{>= IUPACThreshold} will be included in
#' the IUPAC consensus sequence. The default is 0.
#'
#' @section
#' Majority consensus sequence:
#' The most frequently occurring nucleotide.
#' If two or more bases occur with the same frequency,
#' the consensus nucleotide will be randomly selected among these bases.
#'
#' @section
#' IUPAC consensus sequence:
#' The consensus sequence expressed in IUPAC format (i.e. with wobble bases)
#' Note that the IUPAC consensus sequence only
#' takes 'A', 'C', 'G', 'T' and '-' as input. Degenerate bases
#' present in the alignment will be skipped. If a position only contains
#' degenerate/invalid bases, the IUPAC consensus will be \code{NA} at that
#' position.
#'
#' @section Gaps:
#' Proportion of gaps. Gaps are recognized as "-" and ".".
#'
#' @section Identity:
#' Proportion of
#' the most common base. Gaps (-),
#' as well as bases other than A, C, G and T are excluded from the
#' calculation.
#'
#' @section Entropy:
#' Shannon entropy is a measurement of
#' variability. First, for each nucleotide that occurs at a specific position,
#' \code{p*log2(p)}, is calculated, where \code{p} is the proportion of
#' that nucleotide. Then, these values are summarized,
#' followed by multiplication by \code{-1}.
#' A value of \code{0} indicate no variability and a high value
#' indicate high variability.
#' Gaps (-), as well as bases other than
#' A, C, G and T are excluded from the calculation.
#'
#' @return
#' An RprimerProperties object, which contains a tibble (a data frame)
#' with information on majority and IUPAC consensus sequence, gap frequency,
#' nucleotide identity and Shannon entropy.
#'
#' @examples
#' GetAlignmentProperties(example_rprimer_profile)
#' GetAlignmentProperties(example_rprimer_profile, IUPACThreshold = 0.1)
#'
#' @export
getSequenceProperties <- function(x, IUPACThreshold = 0) {
  if (!is.RprimerProfile(x)) {
    stop("'x' must be an RprimerProfile object", call. = FALSE)
  }
  position <- seq_len(ncol(x))
  majority <- majorityConsensus(x)
  IUPAC <- IUPACConsensus(x, threshold = IUPACThreshold)
  gaps <- gapFrequency(x)
  identity <- nucleotideIdentity(x)
  entropy <- shannonEntropy(x)
  sequenceProperties <- tibble::new_tibble( ####################
    list(
      "position" = position,
      "majority" = majority,
      "IUPAC" = IUPAC,
      "gaps" = gaps,
      "identity" = identity,
      "entropy" = entropy
    ),
    nrow = length(position),
    class = "rprimerProperties"
  )
  sequenceProperties
}

# Helpers =====================================================================

#' Majority consensus sequence
#'
#' \code{majorityConsensus} returns the majority consensus sequence of an
#' alignment of DNA sequences.
#'
#' @inheritParams getSequenceProperties
#'
#' @return The majority consensus sequence (a character vector of length n).
#'
#' @keywords internal
#'
#' @noRd
majorityConsensus <- function(x) {
  if (!is.RprimerProfile(x)) {
    stop("'x' must be an RprimerProfile object.", call. = FALSE)
  }
  x <- unclass(x)
  findMostCommonBase <- function(x, y) {
    mostCommon <- rownames(x)[y == max(y)]
    # If there are ties, the most common base will be randomly selected
    if (length(mostCommon > 1)) {
      mostCommon <- sample(mostCommon, 1)
    }
    return(mostCommon)
  }
  # Get the consensus sequence at all positions
  consensus <- apply(x, 2, function(y) findMostCommonBase(x, y))
  consensus <- unname(consensus)
  consensus
}

#' Convert DNA nucleotides into the corresponding IUPAC degenerate base
#'
#' \code{asIUPAC} takes several DNA nucleotides as input,
#' and returns the degenerate base in IUPAC format.
#'
#' @param x A character vector of length one, containing DNA
#' nucleotides (valid bases are A, C, G, T, - and .). Each base must
#' be separated by a comma (,), e.g. 'A,C,G'.
#' Characters other than A, C, G, T, -, ., will be ignored. However, when
#' the input only consist of invalid bases,
#' or if the bases are not separated by ',',
#' \code{asIUPAC} will return NA.
#'
#' @return The corresponding IUPAC base.
#'
#' @examples
#' as_IUPAC("A,C")
#' as_IUPAC("R")
#' as_IUPAC("TG") ## Will return NA since the bases are not separated by comma
#'
#' @keywords internal
#'
#' @noRd
asIUPAC <- function(x) {
  if (!(is.character(x) && length(x) == 1)) {
    stop(
      "'x' must be a character vector of length one, e.g. 'A,C,T'.",
      call. = FALSE
    )
  }
  x <- gsub(" ", "", x)
  x <- unlist(strsplit(x, split = ","), use.names = FALSE)
  x <- x[order(x)]
  x <- unique(x)
  bases <- c("A", "C", "G", "T", "-", ".")
  match <- x %in% bases
  x <- x[!(match == FALSE)] # exclude non accepted nucleotides
  x <- paste(x, collapse = ",")
  IUPAC <- unname(IUPACLookup[x]) # match with lookup table
  IUPAC
}

#' IUPAC consensus sequence
#'
#' \code{IUPACConsensus} returns the IUPAC consensus sequence from an
#' PrprimerProfile object.
#'
#' @inheritParams getSequenceProperties
#'
#' @param threshold
#' Optional. A number (0, 0.2]
#' At each position, all nucleotides with a proportion
#' \code{>= threshold} will be included in
#' the IUPAC consensus sequence. The default is 0.
#'
#' @return The consensus sequence (a character vector of length n).
#'
#' @seealso \code{asIUPAC} for further info on how the IUPAC consensus
#' sequecnce is determined.
#'
#' @keywords internal
#'
#' @noRd
IUPACConsensus <- function(x, threshold = 0) {
  if (!is.RprimerProfile(x)) {
    stop("'x' must be an RprimerProfile object.", call. = FALSE)
  }
  if (!is.double(threshold) || threshold < 0 || threshold > 0.2) {
    stop(paste0(
      "'treshold' must be higher than 0 and less or equal to 0.2. \n
      You've set it to ", threshold, "."
    ))
  }
  x <- unclass(x)
  bases <- c("a", "c", "g", "t", "-", ".")
  x <- x[rownames(x) %in% bases, ]
  basesToInclude <- apply(x, 2, function(y) {
    paste(rownames(x)[y > threshold], collapse = ",")
  })
  basesToInclude <- unname(basesToInclude)
  consensus <- purrr::map_chr(basesToInclude, ~ asIUPAC(.x))
  if (any(is.na(consensus))) {
    warning("The consensus sequence contain NAs. \n
    Try to lower the threshold value.", call. = FALSE)
  }
  consensus
}

#' Gap frequency
#'
#' \code{gapFrequency} returns the gap frequency from an
#' PrprimerProfile object.
#'
#' @inheritParams getSequenceProperties
#'
#' @return The gap frequency (a numeric vector of length n).
#'
#' @keywords internal
#'
#' @noRd
gap_frequency <- function(x) {
  if (!is.RprimerProfile(x)) {
    stop("'x' must be an RprimerProfile object.", call. = FALSE)
  }
  x <- unclass(x)
  if ("-" %in% rownames(x)) {
    gaps <- x[rownames(x) == "-", ] ### also add "."
    gaps <- unname(gaps)
  } else {
    gaps <- rep(0, ncol(x))
  }
  gaps
}

#' Nucleotide identity
#'
#' \code{nucleotideIdentity} returns the nucleotide identity from an
#' PrprimerProfile object.
#'
#' @inheritParams getSequenceProperties
#'
#' @return The nucleotide identity (a numeric vector of length n).
#' The nucleotide identity can range between (0, 1].
#'
#' @keywords internal
#'
#' @noRd
nucleotideIdentity <- function(x) {
  if (!is.RprimerProfile(x)) {
    stop("'x' must be an RprimerProfile object.", call. = FALSE)
  }
  # I want to assess identity based on DNA bases,
  # i.e. ignore gaps and degenerate positions, so I make a subset (s) of x
  # with the rows named a, c, g and t.
  x <- unclass(x)
  bases <- c("A", "C", "G", "T")
  s <- x[rownames(x) %in% bases, ]
  # Calculate relative proportions of the bases in s
  s <- apply(s, 2, function(x) x / sum(x))
  # Find the greatest proportion at each position
  identity <- apply(s, 2, max)
  identity <- unname(identity)
  identity[is.na(identity)] <- 0
  identity
}

#' Shannon entropy
#'
#' \code{shannonEntropy} returns the Shannon entropy from an
#' PrprimerProfile object.
#'
#' @inheritParams getSequenceProperties
#'
#' @return The Shannon entropy (a numeric vector of length n).
#'
#' @keywords internal
#'
#' @noRd
shannonEntropy <- function(x) {
  if (!is.RprimerProfile(x)) {
    stop("'x' must be an RprimerProfile object.", call. = FALSE)
  }
  # I want to assess entrpoy from DNA bases,
  # i.e. ignore gaps and degenerate positions, so I make a subset (s) of x
  # with the rows named a, c, g and t.
  x <- unclass(x)
  bases <- c("A", "C", "G", "T")
  s <- x[rownames(x) %in% bases, ]
  # Calculate proportions of the bases in s
  s <- apply(s, 2, function(x) x / sum(x))
  entropy <- apply(s, 2, function(x) {
    ifelse(x == 0, 0, x * log2(x))
  })
  entropy <- -colSums(entropy)
  entropy <- unname(entropy)
  entropy <- abs(entropy) # abs to avoid -0 (due to neg sums)
  entropy[is.na(entropy)] <- 0
  entropy
}

#' Check if an object is as RprimerProperties-object
#'
#' @param x An 'RprimerProperties'-like object.
#'
#' @return \code{TRUE} or \code{FALSE}.
#'
#' @keywords internal
#'
#' @noRd
is.RprimerProperties <- function(x) inherits(x, "RprimerProperties")
