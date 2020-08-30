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
#' @keywords internal
#'
#' @noRd
nn_split <- function(x) {
  if (!(!is.na(x) && is.character(x) && nchar(x) > 1)) {
    stop("'x' must be a character vector of length one,
      with at least two characters (e.g. 'caaggnt')", call. = FALSE)
  }
  x <- split_sequence(x)
  from <- (seq_along(x) - 1)[-1]
  to <- seq_along(x)[-1]
  nn <- purrr::map2_chr(from, to, function(i, j) {
    paste(x[i:j], collapse = "")
  })
  nn
}

#' Calculate dH or dS of nearest neighbors using lookup tables
#'
#' @param x A matrix (of type character) with nearest-neighbor pairs
#' of DNA sequences (e.g. \code{c('ct', 'tt', 'ta')})
#'
#' @param table The lookup table that should be used.
#' Either 'dH' (entropy) or 'dS' (enthalpy).
#'
#' @details \code{x} cannot contain other characters than
#' 'a', 'c', 'g', 't' and '-'.
#'
#' @return The corresponding values for dH or dS (in cal/M).
#'
#' @examples
#' nn_lookup(c("ac", "cc"), table = "dS")
#'
#' @keywords internal
#'
#' @noRd
nn_lookup <- function(x, table) {
  valid_tables <- c("dH", "dS")
  if (!is.character(table) || !table %in% valid_tables) {
    stop("'table' must be set to either 'dH' or 'dS'", call. = FALSE)
  }
  if (any(grepl("[^acgt]", x))) {
    stop(
      "'x' contains at least one invalid base. \n
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
  result
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
#' @keywords internal
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
#' @param oligos One or more DNA sequences (a character vector).
#'
#' @inheritParams get_oligos
#'
#' @return The melting temperature(s) of x.
#'
#' @examples
#' tm("acggtgcctac")
#' tm(c("acggtgcctac", "acggtggctgc"))
#'
#' @keywords internal
#'
#' @noRd
tm <- function(oligos, conc_oligo = 5e-07, conc_na = 0.05) {
  if (!is.double(conc_oligo) || conc_oligo < 2e-07 || conc_oligo > 2.0e-06) {
    stop("'conc_oligo' must be between
           0.2e-07 M (20 nM) and 2e-06 M (2000 nM)", call. = FALSE)
  }
  if (!is.double(conc_na) || conc_na < 0.01 || conc_na > 1) {
    stop("'conc_na' must be between 0.01 and 1 M", call. = FALSE)
  }
  oligos <- tolower(oligos)
  # Find initiation values
  init_H <- purrr::map_dbl(oligos, ~ init_5end(.x)[["H"]] + init_3end(.x)[["H"]])
  init_S <- purrr::map_dbl(oligos, ~ init_5end(.x)[["S"]] + init_3end(.x)[["S"]])
  nn <- purrr::map(oligos, nn_split)
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
  N <- nchar(oligos[[1]]) - 1 # Number of phosphates
  sumdS <- sumdS + 0.368 * N * log(conc_na)
  tm <- sumdH / (sumdS + gas_constant * log(conc_oligo))
  tm <- tm - 273.15
  tm
}
