#' Get sequence properties
#'
#' \code{sequence_properties} returns sequence information from an alignment
#' of DNA sequences.
#'
#' @param x A nucleotide profile of class 'rprimer_profile'.
#'
#' @param iupac_threshold
#' A number between 0 and 0.2 (the default is 0).
#' At each position, all nucleotides with a proportion
#' higher than or equal to the stated threshold will be included in
#' the iupac consensus sequence.
#'
#' @section
#' Majority consensus sequence:
#' The most frequently occuring nucleotide.
#' If two or more bases occur with the same frequency,
#' the consensus nucleotide will be randomly selected among these bases.
#'
#' @section
#' IUPAC consensus sequence:
#' The consensus sequence expressed in IUPAC format (i.e. with wobble bases)
#' Note that the IUPAC consensus sequence only
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
#' as well as nucleotides other than a, c, g and t, are excluded from the
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
#' a, c, g and t, are excluded from the calculation.
#'
#' @return
#' A tibble (data frame) of class 'rprimer_properties',
#' with information about majority and iupac consensus sequence, gap frequency,
#' nucleotide identity and shannon entropy.
#'
#' \describe{
#'   \item{position}{position in the alignment}
#'   \item{majority}{majority consensus sequence}
#'   \item{iupac}{iupac consensus sequence}
#'   \item{gaps}{proportion of gaps}
#'   \item{identity}{proportion of the most common nucleotide}
#'   \item{entropy}{Shannon entropy}
#' }
#'
#' @examples
#' sequence_properties(example_rprimer_profile)
#'
#' @export
sequence_properties <- function(x, iupac_threshold = 0) {
  if (!inherits(x, "rprimer_profile")) {
    stop("'x' must be an rprimer_profile object.", call. = FALSE)
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
