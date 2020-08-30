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
#' value is 14 and the maximum allowed value is 30.
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
#' allowed value is 16 (default is 4).
#'
#' @param x One or more DNA sequences (a character vector).
#'
#' @param avoid_3end_ta
#' \code{TRUE} or \code{FALSE}.
#' If oligos with t or a at the 3' end
#' should be replaced with \code{NA}
#' (recommended for primers). The default is \code{TRUE}.
#'
#' @param avoid_3end_runs
#' \code{TRUE} or \code{FALSE}.
#' If oligos with more than two runs
#' of the same nucleotide at the 3' end should be replaced with \code{NA}
#' (recommended for primers).
#' The default is \code{TRUE}.
#'
#' @param avoid_gc_rich_3end
#' \code{TRUE} or \code{FALSE}.
#' If oligos with more than three G or C within the last five bases of the 3'
#' end should be replaced with \code{NA}
#' (recommended for primers).
#' The default is \code{TRUE}.
#'
#' @param avoid_5end_g
#' \code{TRUE} or \code{FALSE}.If oligos with g
#' at the 5' end should be replaced with \code{NA}
#' (recommended for probes). The default is \code{FALSE}.
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
#' The melting temperature is calculated using the nearest-neighbour method.
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
#' The tibble contains the following information:
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
#'   \item{tm_majority}{melting temperature (majority sequence),
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
#' avoid_3end_runs = TRUE,
#' avoid_gc_rich_3end = TRUE,
#' avoid_5end_g = FALSE,
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
#' Thermodynamics and NMR of Internal G-T Mismatches in DNA.
#' Biochemistry, 36, 34: 10581–10594
#' (Duplex initiation parameters are from here)
#'
#' SantaLucia, J (1998) A unified view of polymer,
#' dumbell, and oligonucleotide DNA nearest-neighbor thermodynamics.
#' Proc. Natl. Acad. Sci. USA, 95: 1460-1465. (Table values are from here)
#'
#' @export
get_oligos <- function(x,
                       target,
                       length = 18:22,
                       max_gap_frequency = 0.1,
                       max_degenerates = 2,
                       max_degeneracy = 4,
                       avoid_3end_ta = TRUE,
                       avoid_3end_runs = TRUE,
                       avoid_gc_rich_3end = TRUE,
                       avoid_5end_g = FALSE,
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
      avoid_3end_runs = avoid_3end_runs,
      avoid_gc_rich_3end = avoid_gc_rich_3end,
      avoid_5end_g = avoid_5end_g
    )
    oligos
  })
  if (nrow(all_oligos) == 0L)
    stop("No oligos were found.", call. = FALSE)
  # Check match to targets
  all_oligos <- check_match(all_oligos, target)
  all_oligos <- dplyr::arrange(all_oligos, all_oligos$begin)
  all_oligos <- tibble::new_tibble(
    nrow = nrow(all_oligos), class = "rprimer_oligo"
  )
  all_oligos
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
#' get_nmers(c("c", "g", "t", "t", "c", "g"), n = 2)
#'
#' @keywords internal
#'
#' @noRd
get_nmers <- function(x, n = NULL) {
  if (!(is.character(x))) {
    stop("'x' must be a character vector", call. = FALSE)
  }
  if (is.null(n)) {
    n <- round(length(x) / 10)
    if (n == 0) {
      n <- 1
    }
  }
  if (!is.numeric(n) || n < 1 || n > length(x)) {
    stop("'n' must be an integer between 1 and length(x)", call. = FALSE)
  }
  begin <- 1:(length(x) - n + 1)
  end <- begin + n - 1
  nmer <- purrr::map_chr(
    begin, ~ paste(x[begin[[.x]]:end[[.x]]], collapse = "")
  )
  nmer
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
#'
#' @keywords internal
#'
#' @noRd
gc_content <- function(x) {
  if (typeof(x) != "character" || length(x) != 1) {
    stop("'x' must be a character vector of length one", call. = FALSE)
  }
  x <- tolower(x)
  x <- split_sequence(x)
  gc_count <- length(which(x == "c" | x == "g"))
  # Gaps will not be included in the total count
  total_count <- length(which(x == "a" | x == "c" | x == "g" | x == "t"))
  gc <- gc_count / total_count
  gc
}

#' Reverse complement
#'
#' \code{reverse_complement} finds the reverse complement of a DNA seuquence.
#'
#' @param x A DNA sequence (a character vector of length one, e.g. 'cttgg').
#'
#' @details For \code{x}, valid bases are 'a', 'c', 'g', 't', 'r', 'y', 'm',
#' 'k', 's', 'w', n', 'h', 'd', 'v', 'b' and '-'.
#'
#' @return The reverse complement. Non valid bases will return as \code{NA}.
#'
#' @examples
#' reverse_complement("cttgtr")
#'
#' @keywords internal
#'
#' @noRd
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
#' @keywords internal
#'
#' @noRd
running_sum <- function(x, n = NULL) {
  if (!(is.numeric(x))) {
    stop("'x' must be a numeric vector.", call. = FALSE)
  }
  if (is.null(n)) {
    n <- round(length(x) / 10)
    if (n == 0) {
      n <- 1
    }
  }
  if (!is.numeric(n) || n < 1 || n > length(x)) {
    stop("'n' must be a number between 1 and length(x).", call. = FALSE)
  }
  cumul <- c(0, cumsum(x))
  runsum <- cumul[(n + 1):length(cumul)] - cumul[1:(length(cumul) - n)]
  runsum
}

#' Exclude oligos with a specific pattern
#'
#' \code{exclude} replaces oligos with a specific pattern
#' with \code{NA}.
#'
#' @param oligos One or more oligo sequences (a character vector).
#'
#' @param pattern A regular expression.
#'
#' @return A character vector where the oligos with the pattern
#' have been replaced with \code{NA}.
#'
#' @examples
#' exclude(c("cttgttatttt", "cgattctg"), "a$|t$")
#'
#' @keywords internal
#'
#' @noRd
exclude <- function(oligos, pattern) {
  if (!is.character(oligos)) {
    stop("'oligos' must be a character vector.", call. = FALSE)
  }
  if (!is.character(pattern)) {
    stop("'pattern' must be a character vector.", call. = FALSE)
  }
  regex <- pattern
  oligos[grepl(regex, oligos)] <- NA
  oligos
}

#' Exclude oligos
#'
#' \code{exclude oligos} replaces oligos with many consecutive
#' mono- or dinucleotides with \code{NA}. If selected, it also replaces oligos
#' with 5'-end g or 3'-end t or a with \code{NA}.
#'
#' @param oligos One or more oligos (a character vector).
#'
#' @inheritParams get_oligos
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
#'
#' @keywords internal
#'
#' @noRd
exclude_oligos <- function(oligos,
                           avoid_3end_ta = FALSE,
                           avoid_3end_runs = FALSE,
                           avoid_gc_rich_3end = TRUE,
                           avoid_5end_g = FALSE) {
  if (any(!is.logical(
    c(avoid_3end_ta, avoid_3end_runs, avoid_gc_rich_3end, avoid_5end_g)
    ))) {
    stop(
      "'avoid_3end_ta', 'avoid_5end_g', 'avoid_gc_rich_3end' and
      'avoid_3end_runs' must be set to TRUE or FALSE",
      call. = FALSE
    )
  }
  # Remove oligos with at least 4 'runs' of the same dinucleotide ##################
  oligos <- exclude(
    oligos,
    "(at){4,}|(ta){4,}|(ac){4,}|(ca){4,}|(ag){4,}|(ga){4,}|(gt){4,}|(tg){4,}|(cg){4,}|(gc){4,}|(tc){4,}|(ct){4,})"
  )
  # Remove oligos with at least 5 'runs' of the same nucleotide
  oligos <- exclude(oligos, "([a-z])\\1\\1\\1\\1")
  if (avoid_3end_ta == TRUE) {
    # Remove oligos with t or a at the 3' end
    oligos <- exclude(oligos, "a$|t$")
  }
  if (avoid_3end_runs == TRUE) {
    # Remove oligos with at least 3 'runs' of the same nucleotide in 3'
    oligos <- exclude(oligos, "([a-z])\\1\\1$")
  }
  if (avoid_gc_rich_3end == TRUE) {
   # temp <- purrr::map(oligos, ~ split_sequence(.x))
  #  temp <- purrr::map(temp, ~ .x[(length(.x) - 4):length(.x)]) ######################
  #  temp <- purrr::map(temp, ~ paste(.x, collapse = ""))
  #  gc <- purrr::map_dbl(temp, ~ gc_content(.x))
  #  oligos[which(gc < 4/6)] <- NA
  oligos <- oligos
  }
  if (avoid_5end_g == TRUE) {
    # Remove oligos with g at the 5' end
    oligos <- exclude(oligos, "^g")
  }
  oligos
}

#' Exclude unwanted oligos
#'
#' \code{exclude unwanted oligos} replaces oligos with many consecutive
#' mono- or dinucleotides with \code{NA}. If selected, it also replaces oligos
#' with 5'-end g or 3'-end t or a with \code{NA}.
#'
#' @param x A tibble with oligos.
#'
#' @param avoid_3end_ta
#' \code{TRUE} or \code{FALSE}.
#' If oligos with t or a at the 3' end
#' should be replaced with \code{NA}
#' (recommended for primers).
#'
#' @param avoid_3end_runs
#' \code{TRUE} or \code{FALSE}.
#' If oligos with more than two runs
#' of the same nucleotide at the 3' end should be replaced with \code{NA}
#' (recommended for primers).
#'
#' @param avoid_gc_rich_3end
#' \code{TRUE} or \code{FALSE}.
#' If oligos with more than three G or C within the last five bases of the 3' end
#' should be replaced with \code{NA}
#' (recommended for primers).
#'
#' @param avoid_5end_g
#' \code{TRUE} or \code{FALSE}.If oligos with g
#' at the 5' end should be replaced with \code{NA}
#'
#' @details An oligo is replaced with \code{NA} if
#' - It has more than three runs of the same dinucleotide (e.g. 'tatatata')
#' - It has more than four runs of the same nucleotide
#'
#' @return A tibble where unwanted oligos
#' has been replaced with \code{NA}.
#'
#' @keywords internal
#'
#' @noRd
exclude_unwanted_oligos <- function(x,
                                    avoid_3end_ta,
                                    avoid_3end_runs,
                                    avoid_gc_rich_3end,
                                    avoid_5end_g) {
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
  x
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
#'
#' @keywords internal
#'
#' @noRd
count_degenerates <- function(x) {
  if (typeof(x) != "character") {
    stop("'x' must be a character vector.", call. = FALSE)
  }
  if (grepl("[^acgtrymkswnhdvb-]", x)) {
    stop("'x' contains at least one invalid base. \n
      Valid bases are 'a', 'c', 'g', 't', 'r', 'y', 'm', 'k', 's', 'w',
      'n', 'h', 'd', 'v', 'b' and '-'",
         call. = FALSE
    )
  }
  nt <- c("a", "c", "g", "t", "-")
  x <- split_sequence(x)
  count <- length(x[!x %in% nt])
  count
}

#' Count the degeneracy of a DNA sequence
#'
#' \code{count_degenerates} counts the number of unique sequences of
#' a DNA sequence with degenerate bases.
#'
#' @param x a DNA sequence (a character vector of length one, e.g. 'cttgg').
#'
#' @details
#' Valid bases for \code{x} are 'a', 'c', 'g', 't', 'r', 'y', 'm',
#' 'k', 's', 'w', n', 'h', 'd', 'v', 'b' and '-'.
#'
#' @return The number of unique sequences of x (an integer).
#'
#' @examples
#' count_degeneracy("cttnra")
#'
#' @keywords internal
#'
#' @noRd
count_degeneracy <- function(x) {
  if (typeof(x) != "character") {
    stop("'x' must be a character vector.", call. = FALSE)
  }
  if (grepl("[^acgtrymkswnhdvb-]", x)) {
    stop("'x' contains at least one invalid base. \n
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
  degeneracy
}

#' Generate oligos of a specific length
#'
#' @inheritParams get_oligos
#'
#' @return A tibble with all possible oligos.
#'
#' @keywords internal
#'
#' @noRd
generate_oligos <- function(x,
                            oligo_length = 20,
                            max_gap_frequency = 0.1,
                            max_degenerates = 2,
                            max_degeneracy = 4) {
  if (!is.rprimer_properties(x)) {
    stop(
      "'x' must be an rprimer_properties object.",
      call. = FALSE
    )
  }
  if (!(min(oligo_length) >= 14 && max(oligo_length) <= 30)) {
    stop("'oligo_length' must be between 14 and 30", call. = FALSE)
  }
  if (!(max_gap_frequency >= 0 && max_gap_frequency <= 1)) {
    stop("'max_gap_frequency' must be between 0 and 1", call. = FALSE)
  }
  if (!(max_degenerates <= 6 && max_degenerates >= 0)) {
    stop("'max_degenerates' must be between 0 and 6", call. = FALSE)
  }
  if (!(max_degeneracy >= 1 && max_degeneracy <= 16)) {
    stop("'max_degeneracy' must be between 1 and 16", call. = FALSE)
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
  gap_bin <- ifelse(x$gaps > max_gap_frequency, 1L, 0L) ###################
  gap_penalty <- running_sum(gap_bin, n = oligo_length)
  oligos <- tibble::tibble(
    begin, end, length, majority, iupac,
    majority_rc, iupac_rc, degenerates, degeneracy, gap_penalty
  )
  # Identify and exclude oligos that are duplicated
  unique_oligos <- match(oligos$majority, unique(oligos$majority))
  oligos <- oligos[unique_oligos, ]
  # Exclude oligos with too high gap frequency
  oligos <- oligos[oligos$gap_penalty == 0, ]
  # Exclude oligos with too many degenerate bases
  oligos <- oligos[oligos$degenerates <= max_degenerates, ]
  # Exclude oligos with too high degeneracy
  oligos <- oligos[oligos$degeneracy <= max_degeneracy, ]
  # Remove the gap_penalty column
  oligos <- dplyr::select(oligos, -gap_penalty)
  oligos
}

#' Calculate GC content and tm of oligos
#'
#' @param oligos A tibble with oligos.
#'
#' @inheritParams get_oligos
#'
#' @return A tibble (a data frame) with oligo candidates.
#'
#' @keywords internal
#'
#' @noRd
add_gc_tm <- function(oligos,
                      gc_range = c(0.45, 0.55),
                      tm_range = c(48, 70),
                      conc_oligo = 5e-07,
                      conc_na = 0.05) {
  if (!(min(gc_range) >= 0 && max(gc_range) <= 1)) {
    stop(
      "'gc_range' must be between 0 and 1, e.g. c(0.45, 0.65)", call. = FALSE
    )
  }
  if (!(min(tm_range) >= 20 && max(tm_range) <= 90)) {
    stop(
      "'tm_range' must be between 20 and 90, e.g. c(55, 60)", call. = FALSE
    )
  }
  # Calculate GC content of all majority oligos
  gc_majority <- purrr::map_dbl(oligos$majority, ~ gc_content(.x))
  oligos <- tibble::add_column(oligos, gc_majority)
  # Exclude oligos with GC content outside the stated thresholds
  oligos <- oligos[oligos$gc_majority >= min(gc_range), ]
  oligos <- oligos[oligos$gc_majority <= max(gc_range), ]
  # Calculate Tm of all majority oligos
  if (nrow(oligos) > 0L) {
    tm_majority <- tm(oligos$majority, conc_oligo = conc_oligo, conc_na = conc_na)
  } else {
    tm_majority <- NA
  }
  oligos <- tibble::add_column(oligos, tm_majority)
  # Exclude oligos with Tm outside the stated thresholds
  oligos <- oligos[oligos$tm_majority >= min(tm_range), ]
  oligos <- oligos[oligos$tm_majority <= max(tm_range), ]
  oligos
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
#'
#' @keywords internal
#'
#' @noRd
make_regex <- function(x) {
  if (typeof(x) != "character" || length(x) != 1) {
    stop("'x' must be a character vector of length one", call. = FALSE)
  }
  if (grepl("[^acgtryswkmbdhvn-]", x)) {
    stop("'x' contains at least one invalid base. \n
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
  regx
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
#' @keywords internal
#'
#' @noRd
check_match <- function(x, y) {
  if (!is.rprimer_alignment(y)) {
    stop("'y' must be an rprimer_alignment object.", call. = FALSE)
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
  x
}

#' Check if an object has class attribute rprimer_oligo
#'
#' @param x An rprimer_oligo-like object.
#'
#' @return \code{TRUE} or \code{FALSE}.
#'
#' @keywords internal
#'
#' @noRd
is.rprimer_oligo <- function(x) inherits(x, "rprimer_oligo")
