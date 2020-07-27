get_oligos <- function(x, length = 18:22, max_gap_frequency = 0.1, max_degenerates = 2, max_degeneracy = 4, avoid_3end_ta = TRUE,
                       avoid_5end_g = FALSE, avoid_3end_runs = TRUE, gc_range = c(0.45, 0.55), tm_range = c(48, 70), conc_oligo = 5e-07, conc_na = 0.05) {





  #generate_oligos
  #exclude_unwanted_oligos
  #7check_match

  if (nrow(all_oligos) == 0L)
    stop("No oligos were found.", call. = FALSE)
  all_oligos <- dplyr::arrange(all_oligos, begin)
  all_oligos <- tibble::new_tibble(all_oligos, nrow = nrow(all_oligos), class = "rprimer_oligo")
  return(all_oligos)
}

#' Generate oligos of a specific length
#'
#' @param x An object of class 'rprimer_sequence_properties'
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
  max_degeneracy = 4,
  gc_range =c(0.45, 0.55),
  tm_range = c(48, 70),
  conc_oligo = 5e-07,
  conc_na = 0.05
  ) {
  if (!inherits(x, "rprimer_sequence_properties")) {
    stop(
      "An rprimer_sequence_properties object is expected for x.", call. = FALSE
    )
  }
  if (!(min(oligo_length) >= 6 && max(oligo:length) <= 30)) {
    stop("oligo_length must be between 6 and 30", call. = FALSE)
  }
  if (!(max_gap_frequency >= 0 && max_gap_frequency <= 1)) {
    stop("max_gap_frequency must be between 0 and 1", call. = FALSE)
  }
  if (!(max_degenerates <= 30 && max_degenerates >= 0)) {
    stop("max_degenerates must be between 0 and 6", call. = FALSE)
  }
  if (!(max_degeneracy >= 1 && max_degeneracy <= 64)) {
    stop("max_degeneracy must be between 1 and 64", call. = FALSE)
  }
  if (!(min(gc_range) >= 0 && max(gc_range) <= 1)) {
    stop("gc_range must be between 0 and 1, e.g. c(0.45, 0.65)", call. = FALSE)
  }
  if (!(min(tm_range) >= 20 && max(tm_range) <= 90)) {
    stop("tm_range must be between 20 and 90, e.g. c(55, 60)", call. = FALSE)
  }

  # Find all possible oligos of length y
  majority <- get_nmers(x$majority, n = oligo_length)
  iupac <- get_nmers(x$iupac, n = oligo_length)
  majority_rc <- purrr::map_chr(majority, ~reverse_complement(.x))
  iupac_rc <- purrr::map_chr(iupac, ~reverse_complement(.x))
  degenerates <- purrr::map_int(iupac, ~count_degenerates(.x))
  degeneracy <- purrr::map_dbl(iupac, ~count_degeneracy(.x))
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
  # Calculate GC content of all majority oligos
  gc_majority <- purrr::map_dbl(oligos$majority, ~gc_content(.x))
  oligos <- tibble::add_column(oligos, gc_majority)
  # Exclude oligos with GC content outside the stated thresholds
  oligos <- oligos[oligos$gc_majority >= min(gc_range), ]
  oligos <- oligos[oligos$gc_majority <= max(gc_range), ]
# Calculate Tm of all majority oligos
  tm_majority <- tm(oligos$majority, conc_oligo = conc_oligo, conc_na = conc_na)
  oligos <- tibble::add_column(oligos, tm_majority)
  # Exclude oligos with Tm outside the stated thresholds
  oligos <- oligos[tm_majority >= min(tm_range), ]
  oligos <- oligos[tm_majority <= max(tm_range), ]
  return(oligos)
}

#filter_oligos <- function(oligos) {
#  }



all_oligos <- purrr::map_dfr(length, function(y) {
  # Find all possible oligos of length y


  # Calculate GC content of all majority oligos
  gc_majority <- purrr::map_dbl(oligos$majority, ~gc_content(.x))
  oligos <- tibble::add_column(oligos, gc_majority)
  # Exclude oligos with GC content outside the stated thresholds
  oligos <- oligos[which(gc_majority >= min(gc_range) & gc_majority <= max(gc_range)), ]
  # Calculate Tm of all majority oligos
  tm_majority <- tm(oligos$majority)
  oligos <- tibble::add_column(oligos, tm_majority)
  # Exclude oligos with Tm outside the stated thresholds
  oligos <- oligos[which(tm_majority >= min(tm_range) & tm_majority <= max(tm_range)), ]
  # Indentify non-complex oligos and (if specified) oligos with 3' t-a or 5' g, and set to NA
  oligos$majority <- exclude_unwanted_oligos(oligos$majority, avoid_3end_ta = avoid_3end_ta, avoid_5end_g = avoid_5end_g,
                                             avoid_3end_runs = avoid_3end_runs)
  oligos$majority_rc <- exclude_unwanted_oligos(oligos$majority_rc, avoid_3end_ta = avoid_3end_ta, avoid_5end_g = avoid_5end_g,
                                                avoid_3end_runs = avoid_3end_runs)
  oligos$iupac[which(is.na(oligos$majority))] <- NA
  oligos$iupac_rc[which(is.na(oligos$majority_rc))] <- NA
  # Identify oligos where both the sense and antisense sequence is NA
  invalid_oligos <- dplyr::intersect(which(is.na(oligos$majority)), which(is.na(oligos$majority_rc)))
  if (length(invalid_oligos) > 0L)
    oligos <- oligos[-invalid_oligos, ]
  return(oligos)
})
