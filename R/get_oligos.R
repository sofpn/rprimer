#' Get oligos from sequence properties
#'
#' \code{get_oligos} identifies oligos (primers and probes) from
#' sequence properties.
#'
#' @param x an object of class 'rprimer_properties'.
#'
#' @param target
#' The intended target. An alignment of DNA sequences
#' (an object of class 'rprimer_alignment').
#'
#' @param max_gap_frequency
#' Maximum allowed gap frequency.
#' A number between 0 and 1 (default is 0.1, which means that
#' positions with a gap frequency of maximum 0.1 will be
#' considered as an oligo region).
#'
#' @param length
#' Oligo length. The minimum allowed
#' value is 6 and the maximum allowed value is 30.
#' The default is \code{18:22}.
#'
#' @param max_degenerates
#' Maximum number of degenerate positions.
#' The minimum allowed value is 0 and the maximum
#' allowed value is 6 (default is 2).
#'
#' @param max_degeneracy
#' Maximum number of degenerate variants.
#' The minimum allowed value is 1 and the maximum
#' allowed value is 64 (default is 4).
#'
#' @param x One or more DNA sequences (a character vector).
#'
#' @param avoid_3end_ta
#' \code{TRUE} or \code{FALSE}.
#' If oligos with t or a at the 3' end
#' should be replaced with \code{NA}. The default is \code{FALSE}.
#'
#' @param avoid_5end_g
#' \code{TRUE} or \code{FALSE}.If oligos with g
#' at the 5' end should be replaced with \code{NA}. The default is \code{FALSE}.
#'
#' @param avoid_3end_runs
#' \code{TRUE} or \code{FALSE}.
#' If oligos with more than two runs
#' of the same nucleotide at the 3' end should be replaced with \code{NA}.
#' The default is \code{FALSE}.
#'
#' @param gc_range
#' Oligo GC-content range (proportion) Can range between
#' 0 and 1. The default is \code{c(0.45, 0.55)}.
#'
#' @param tm_range
#' The Tm-range of each oligo. Can range between 20 and 90.
#' The default is \code{c(48, 70)}.
#'
#' @param conc_oligo
#' The concentration of oligonucleotide in M,
#' ranging from 0.2e-07 M (20 nM) to 2e-06 M (2000 nM).
#' The default value is 5e-07 M (500 nM) (for Tm calculation)
#'
#' @param conc_na
#' The sodium ion concentration in M, ranging
#' from 0.01 M to 1 M. The default value is 0.05 M (50 mM)
#' (for Tm calculation).
#'
#' @section Excluded oligos:
#' The function excludes oligos with:
#' * More than than three consecutive runs of the same dinucleotide
#' (e.g. 'tatatata')
#' * More than four consecutive runs of the
#' same nucleotide
#' * Oligos that are duplicated
#'
#' @section Tm:
#' The melting temperature is calculated using the nearest-neigbour method.
#'
#' * Oligos are not expected to be self-complementary, so no symmetry
#' correction is done
#' * The oligo concentration is assumed to be much higher
#' than the target concentration
#'
#' See references for table values and equations.
#'
#' @section Note:
#' GC-content and Tm are calculated based on the majority oligos, and
#' may thus be misleading for degenerate (iupac) oligos.
#'
#' @return
#' A tibble (a data frame) of class 'rprimer_oligo', with all oligo
#' candidates. An error message will return if no oligos are found.
#'
#' \describe{
#'   \item{begin}{position where the oligo begins}
#'   \item{end}{position where the oligo ends}
#'   \item{length}{length of the oligo}
#'   \item{majority}{majority sequence}
#'   \item{iupac}{iupac sequence (i.e. with degenerate bases)}
#'   \item{majority_rc}{majority sequence, reverse complement}
#'   \item{iupac_rc}{iupac sequence, reverse complement}
#'   \item{degenerates}{number of degenerate bases}
#'   \item{degeneracy}{number of variants}
#'   \item{gc_majority}{gc-content (majority sequence), proportion}
#'   \item{tm_majority}{melting temperature (majority sequcence),
#'   degrees Celcius}
#'   \item{pm_majority}{proportion of perfectly matching sequences,
#'   (to the target alignment) majority sequence}
#'   \item{pm_iupac}{proportion of perfectly matching sequences,
#'   iupac sequence}
#'   \item{match_matrix}{a logical matrix describing which sequences
#'   in the target alignment the oligo matches perfectly to}
#' }
#'
#' @examples
#' get_oligos <- function(
#' example_rprimer_properties,
#' target_alignment = example_rprimer_alignment,
#' length = 18:22,
#' max_gap_frequency = 0.1,
#' max_degenerates = 2,
#' max_degeneracy = 4,
#' avoid_3end_ta = TRUE,
#' avoid_5end_g = FALSE,
#' avoid_3end_runs = TRUE,
#' gc_range = c(0.45, 0.55),
#' tm_range = c(48, 70)
#' )
#'
#' @references
#' Tm-calculation:
#'
#' SantaLucia, J, et al. (1996)
#' Improved Nearest-Neighbor Parameters for Predicting DNA Duplex Stability.
#' Biochemistry, 35: 3555-3562 (Formula and salt correction method
#' are from here)
#'
#' Allawi, H. & SantaLucia, J. (1997)
#' Thermodynamics and NMR of Internal G·T Mismatches in DNA.
#' Biochemistry, 34: 10581?\200?10594 ###########################################################
#' (Duplex initiation parameters are from here)
#'
#' SantaLucia, J (1998) A unified view of polymer,
#' dumbell, and oligonucleotide DNA nearest-neighbor thermodynamics.
#' Proc. Natl. Acad. Sci. USA, 95: 1460-1465. (Table values are from here)
#'
#' @export
get_oligos <- function(
  x,
  target,
  length = 18:22,
  max_gap_frequency = 0.1,
  max_degenerates = 2,
  max_degeneracy = 4,
  avoid_3end_ta = TRUE,
  avoid_5end_g = FALSE,
  avoid_3end_runs = TRUE,
  gc_range = c(0.45, 0.55),
  tm_range = c(48, 70),
  conc_oligo = 5e-07,
  conc_na = 0.05
) {
  all_oligos <- purrr::map_dfr(length, function(y) {
    oligos <- generate_oligos(
      x,
      oligo_length = y,
      max_gap_frequency = max_gap_frequency,
      max_degenerates = max_degenerates,
      max_degeneracy = max_degeneracy
    )
    oligos <- add_gc_tm(
      oligos,
      gc_range = gc_range,
      tm_range = tm_range,
      conc_oligo = conc_oligo,
      conc_na = conc_na
    )
    oligos <- exclude_unwanted_oligos(
      oligos,
      avoid_3end_ta = avoid_3end_ta,
      avoid_5end_g = avoid_5end_g,
      avoid_3end_runs = avoid_3end_runs
    )
    oligos
  })
  if (nrow(all_oligos) == 0L)
    stop("No oligos were found.", call. = FALSE)
  # Check match to targets
  all_oligos <- check_match(all_oligos, target)
  all_oligos <- dplyr::arrange(all_oligos, all_oligos$begin)
  all_oligos <- tibble::new_tibble(
    all_oligos, nrow = nrow(all_oligos), class = "rprimer_oligo"
  )
  all_oligos
}
