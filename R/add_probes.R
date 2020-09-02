#' Add probes to (RT)-PCR assays
#'
#' \code{add_probes} adds probes to (RT)-PCR assays.
#'
#' @param x Assays to add probes to (an object of class 'rprimer_assay').
#'
#' @param y Candidate probes (an object of class 'rprimer_oligo').
#'
#' @param tm_difference
#' Acceptable Tm difference between primers and probe.
#' The minimum allowed value is -20 and the maximum allowed value is 20.
#' The default are \code{c(0, 20)}.
#'
#' @details
#' The Tm-difference is calculated by subtracting the
#' Tm of the probe with the average Tm of the
#' primer pair. Thus, a negative Tm-difference
#' means that the Tm of the probe is lower than the average Tm of the
#' primer pair.
#'
#' Note that the Tm-difference is calculated from the majority oligos, and
#' may thus be misleading for degenerate (iupac) oligos.
#'
#' @return Assays with probes. A tibble (a data frame) of
#' class 'rprimer_assay'.
#'
#' ##### Describe more here #####
#'
#' @examples
#' add_probes(
#' example_rprimer_assay,
#' example_rprimer_oligo,
#' tm_difference = c(-2, 10)
#' )
#'
#' @seealso get_assays
#'
#' @export
add_probes <- function(x, y, tm_difference = c(0, 20)) {
  if (!is.rprimer_assay(x)) {
    stop("'x' must be an rprimer_assay object.", call. = FALSE)
  }
  if (any(grepl("_pr$", names(x)))) {
    stop("'x' appears to have probes already.", call. = FALSE)
  }
  if (!inherits(y, "rprimer_oligo")) {
    stop("'y' must be an rprimer_oligo object.", call. = FALSE)
  }
  if (!(min(tm_difference) >= -20 && max(tm_difference) <= 20)) {
    stop(
      "'tm_difference' must be between -20 and 20, e.g. c(-1, 5)",
      call. = FALSE
    )
  }
  assays <- x
  probes <- y
  # For each assay, take all probes that bind within the assay region
  probe_candidates <- purrr::map(seq_len(nrow(assays)), function(i) {
    # The probe has to begin after the fwd primer ends
    # and we want at least one base between (hence the + 2)
    from <- assays$end_fwd[[i]] + 2
    # The probe has to end before the rev primer begins
    # and we want at least one base in-between
    to <- assays$begin_rev[[i]] - 2
    # Take all probes that begin and end 'within' the assay region
    probe <- probes[which(probes$begin >= from & probes$end <= to), ]
    probe
  })
  # For each assay, we check how many probe candidates we have
  number_of_probes <- purrr::map_int(probe_candidates, nrow)
  # Then, we need to pick all the assays that can harbor a probe.
  # For assays that has n number of probes, we need to
  # repeat that row n times.
  rows_to_select <- purrr::map(
    seq_along(number_of_probes), ~rep(.x, number_of_probes[[.x]])
  )
  rows_to_select <- unlist(rows_to_select, use.names = FALSE)
  assays <- assays[rows_to_select, ]
  if (nrow(assays) == 0L)
    stop("No assays could be generated.", call. = FALSE)
  # Now, we can make a data frame of the probe candidates
  probe_candidates <- do.call("rbind", probe_candidates)
  # Check the sense of the probe
  # if the plus sense probe is ok, (i.e. not NA), we take that one.
  sense <- ifelse(!is.na(probe_candidates$majority), "pos", "neg")
  probe_candidates$majority <- ifelse(
    sense == "pos", probe_candidates$majority, probe_candidates$majority_rc
  )
  probe_candidates$iupac <- ifelse(
    sense == "pos", probe_candidates$iupac, probe_candidates$iupac_rc
  )
  probe_candidates <- tibble::add_column(probe_candidates, sense)
  drop <- c("majority_rc", "iupac_rc")
  probe_candidates <- probe_candidates[, !names(probe_candidates) %in% drop]
  names(probe_candidates) <- paste0(names(probe_candidates), "_pr")
  assays <- dplyr::bind_cols(assays, probe_candidates)
  assays$total_degeneracy <- assays$total_degeneracy + probe_candidates$degeneracy_pr
  tm_difference_primer_probe <- purrr::map_dbl(
    seq_len(nrow(assays)), function(x) {
      assays$tm_majority_pr[[x]] - mean(
        assays$tm_majority_fwd[[x]], assays$tm_majority_rev[[x]]
      )
    })
  assays <- tibble::add_column(
    assays, tm_difference_primer_probe, .after = "tm_difference_primer"
  )
  assays <- assays[assays$tm_difference_primer_probe >= min(tm_difference), ]
  assays <- assays[assays$tm_difference_primer_probe <= max(tm_difference), ]
  if (nrow(assays) == 0L)
    stop("No assays could be generated.", call. = FALSE)
  assays <- add_probe_to_match_matrix(assays)
  assays <- tibble::new_tibble(assays, class = "rprimer_assay")
  assays
}

# Helpers =====================================================================

#' Add probes to match matrices
#'
#' @param x A tibble with assays (with probes).
#'
#' @return A tibble with assays, with information on perfect matches
#' for forward and reverse primers and probes.
#'
#' @keywords internal
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
    if (!is.matrix(match)) {
      match <- setNames(match) #########################################
      match <- as.matrix(match)
    }
    match <- match[, -(5:6)]
    if (!is.matrix(match)) {
      match <- as.matrix(match)
    }
    majority_all_matr <- match[, c(1, 3, 5)]
    if (!is.matrix(majority_all_matr)) {
      majority_all_matr <- as.matrix(majority_all_matr)
    }
    majority_all <- rowSums(majority_all_matr)
    iupac_all_matr <- match[, c(2, 4, 6)]
    if (!is.matrix(iupac_all_matr)) {
      iupac_all_matr <- as.matrix(iupac_all_matr)
    }
    iupac_all <- rowSums(iupac_all_matr)
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
  x
}
