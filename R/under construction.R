
#' Check if oligos or assays matches their targets (generic)
#'
#' \code{check_match} checks if oligos or assays matches with their
#' intended target sequences.
#'
#' @param x An object of class
#' 'rprimer_oligo' or 'rprimer_assay'.
#'
#' @param y The intended target. An alignment of DNA sequences
#' (an object of class 'rprimer_alignment').
#'
#' @return A tibble (
#' a data frame) of class 'rprimer_oligo' or 'rprimer_assay',
#' with columns describing the proportion of perfectly matching
#' sequences for each oligo or assay, and a column named
#' 'match_report', which contains a matrix with information about
#' which sequences the oligo/assay matches perfectly to.
#'
#' @examples
#' check_match(example_rprimer_oligo, example_rprimer_alignment)
#' check_match(example_rprimer_assay, example_rprimer_alignment)
#'
#' @export
check_match <- function(x, y) {
  if (!inherits(y, "rprimer_alignment")) {
    stop("An rprimer_alignment object is expected for y.", call. = FALSE)
  }
  UseMethod("check_match")
}

#' @describeIn check_match Check match of an object of class
#' 'rprimer_oligo'.
#' @export
check_match.rprimer_oligo <- function(x, y) {
  # If the positive-sense sequence is NA, then take the negative sense and make a reverse complement Make regular
  # expressions of majority and iupac sequences
  majority <- ifelse(!is.na(x$majority), x$majority, purrr::map_chr(x$majority_rc, reverse_complement))
  majority <- purrr::map_chr(majority, make_regex)
  iupac <- ifelse(!is.na(x$iupac), x$iupac, purrr::map_chr(x$iupac_rc, reverse_complement))
  iupac <- purrr::map_chr(iupac, make_regex)
  # Shorten the sequence names to only accession numbers
  names(y) <- purrr::map_chr(names(y), truncate_name)
  # Make a matrix that describes exactly which sequences the oligo matches perfectly to
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
  x <- tibble::new_tibble(x, nrow = nrow(x), class = "rprimer_oligo")
  return(x)
}

#' @describeIn check_match Check match of an object of class
#' 'rprimer_assay'.
#' @export
check_match.rprimer_assay <- function(x, y) {
  # Shorten the names of y to only accession numbers
  names(y) <- purrr::map_chr(names(y), truncate_name)
  # Get match matrices for each assay
  primer_match_matrix <- check_primer_match(x, y)
  if (any(grepl("_pr$", names(x)))) {
    probe_match_matrix <- check_probe_match(x, y)
    match_matrix <- purrr::map(seq_len(nrow(x)), function(x) {
      match <- cbind(primer_match_matrix[[x]], probe_match_matrix[[x]])
      return(match)
    })
  } else match_matrix <- primer_match_matrix
  match_matrix <- purrr::map(match_matrix, function(x) {
    majority_cols <- seq(1, ncol(x), 2)
    iupac_cols <- seq(2, ncol(x), 2)
    majority_all <- rowSums(x[, majority_cols])
    iupac_all <- rowSums(x[, iupac_cols])
    majority_all <- ifelse(majority_all == length(majority_cols), TRUE, FALSE)
    iupac_all <- ifelse(iupac_all == length(iupac_cols), TRUE, FALSE)
    x <- cbind(x, majority_all, iupac_all)
    return(x)
  })
  names(match_matrix) <- paste0(x$majority_fwd, "_", x$majority_rev, "_", x$majority_pr)
  match_percentage <- purrr::map(match_matrix, colMeans)
  match_percentage <- do.call("rbind", match_percentage)
  match_percentage <- tibble::as_tibble(match_percentage)
  names(match_percentage) <- paste0("pm_", names(match_percentage))
  if (any(grepl("^pm_", names(x)))) {
    stop("matches have already been checked in x", call. = FALSE)
  }
  x <- dplyr::bind_cols(x, match_percentage)
  match_matrix <- tibble::tibble(match_matrix)
  x <- dplyr::bind_cols(x, match_matrix)
  x <- tibble::new_tibble(x, nrow = nrow(x), class = "rprimer_assay")
  return(x)
}
