#' Get (RT)-PCR assays from oligos
#'
#' \code{get_assays} combines forward and reverse primers to (RT)-PCR assays.
#'
#' @param x An object of class 'rprimer_oligo'.
#'
#' @param length
#' Amplicon length. Can range from 40 to 5000 base pairs. The
#' default is \code{65:120}.
#'
#' @param max_tm_difference
#' Maximum Tm difference (in C) between the two primers
#' (absolute value). A number between 0 and 30. The default is 1.
#'
#' @details
#' The Tm-difference is calculated from the majority oligos, and
#' may thus be misleading for degenerate (iupac) oligos.
#'
#' @return
#' A tibble (a data frame) of class 'rprimer_assay' with all candidate
#' assays. An error message will return if no assays are found.
#'
#'  \describe{
#'   \item{begin}{position where the assay begins}
#'   \item{end}{position where the assay ends}
#'   \item{amplicon_length}{length of the amplicon}
#'   \item{tm_difference_primer}{difference in melting temperature between
#'   the forward and reverse primer, absolute value, degrees Celcius}
#'   \item{total_degeneracy}{total number of oligos in the assay}
#'   \item{begin_fwd}{position where the forward primer begins}
#'   \item{end_fwd}{position where the reverse primer ends}
#'   \item{length_fwd}{length of the forward primer}
#'   \item{majority_fwd}{majority sequence of the forward primer}
#'   \item{iupac_fwd}{iupac sequence (i.e. with degenerate bases)
#'   of the forward primer}
#'   \item{degenerates_fwd}{number of degenerate bases of the forward primer}
#'   \item{degeneracy_fwd}{number of variants of the forward primer}
#'   \item{gc_majority_fwd}{gc-content of the forward primer
#'   (majority sequence), proportion}
#'   \item{tm_majority_fwd}{melting temperature of the forward primer
#'   (majority sequcence), degrees Celcius}
#'   \item{pm_majority_fwd}{proportion of sequences in the target alignment
#'   that matches perfectly with the forward primer, majority sequence}
#'   \item{pm_iupac_fwd}{proportion of sequences in the target alignment
#'   that matches perfectly with the forward primer, iupac sequence}
#'   \item{begin_rev}{position where the reverse primer begins}
#'   \item{end_rev}{position where the reverse primer ends}
#'   \item{length_rev}{length of the reverse primer}
#'   \item{majority_rev}{majority sequence of the reverse primer}
#'   \item{iupac_rev}{iupac sequence (i.e. with degenerate bases)
#'   of the reverse primer}
#'   \item{degenerates_rev}{number of degenerate bases of the reverse primer}
#'   \item{degeneracy_rev}{number of variants of the reverse primer}
#'   \item{gc_majority_rev}{gc-content of the reverse primer
#'   (majority sequence), proportion}
#'   \item{tm_majority_rev}{melting temperature of the reverse primer
#'   (majority sequcence), degrees Celcius}
#'   \item{pm_majority_rev}{proportion of sequences in the target alignment
#'   that matches perfectly with the reverse primer, majority sequence}
#'   \item{pm_iupac_rev}{proportion of sequences in the target alignment
#'   that matches perfectly with the reverse primer, iupac sequence}
#'   \item{pm_majority_all}{proportion of sequences in the target alignment
#'   that matches perfectly with both the forward and reverse primer,
#'   majority sequence}
#'   \item{pm_iupac_all}{proportion of sequences in the target alignment
#'   that matches perfectly with both the forward and reverse primer,
#'   iupac sequence}
#'   \item{match_matrix}{a logical matrix describing which sequences
#'   in the target alignment the assay matches perfectly to}
#' }
#'
#' @examples
#' get_assays(example_rprimer_oligo, length = 60:150, max_tm_difference = 1.5)
#'
#' @export
get_assays <- function(x, length = 65:120, max_tm_difference = 1) {
  if (!is.rprimer_oligo(x)) {
    stop("'x' must be an rprimer_oligo object.", call. = FALSE) 
  }
  if (!(max_tm_difference > 0 && max_tm_difference < 30)) {
    stop("'max_tm_difference' must be between 0 and 30", call. = FALSE)
  }
  if (!(min(length) >= 40 && max(length) <= 5000)) {
    stop("'length' must be between 40 and 5000", call. = FALSE)
  }
  # Get all pontential candidates for fwd and rev primers
  fwd <- x[!is.na(x$majority), ]
  rev <- x[!is.na(x$majority_rc), ]
  # Get all possible combinations of fwd and rev primers
  combinations <- expand.grid(
    fwd$majority, rev$majority_rc, stringsAsFactors = FALSE
  )
  names(combinations) <- c("fwd_majority", "rev_majority")
  # Get indexes of the fwd and rev primer sequences
  index_fwd <- match(combinations$fwd_majority, x$majority)
  index_rev <- match(combinations$rev_majority, x$majority_rc)
  # Make two datasets of x, one for fwd and one for rev
  fwd <- x[index_fwd, ]
  rev <- x[index_rev, ]
  # Add a tag on colnames before combining the two datasets
  colnames(fwd) <- paste0(colnames(fwd), "_fwd")
  colnames(rev) <- paste0(colnames(rev), "_rev")
  # Combine the two datasets
  assays <- dplyr::bind_cols(fwd, rev)
  assays <- tibble::as_tibble(assays)
  # Add amplicon length, tm difference and total degeneracy to the data
  amplicon_length <- assays$end_rev - assays$begin_fwd + 1
  amplicon_length <- as.integer(amplicon_length)
  tm_difference_primer <- assays$tm_majority_fwd - assays$tm_majority_rev
  tm_difference_primer <- abs(tm_difference_primer)
  begin <- assays$begin_fwd
  end <- assays$end_rev
  total_degeneracy <- assays$degeneracy_fwd + assays$degeneracy_rev
  assays <- tibble::add_column(
    assays, begin, end, amplicon_length,
    tm_difference_primer, total_degeneracy, .before = "begin_fwd"
  )
  # Drop columns that we do no longer need
  drop <- c("majority_rc_fwd", "iupac_rc_fwd", "majority_rev", "iupac_rev")
  assays <- assays[, !(names(assays) %in% drop)]
  # Rename columns
  names(assays)[grep("_rc", names(assays))] <- c("majority_rev", "iupac_rev")
  # Collect assays with desired amplicon length and tm difference
  assays <- assays[assays$amplicon_length >= min(length), ]
  assays <- assays[assays$amplicon_length <= max(length), ]
  assays <- assays[assays$tm_difference_primer <= max_tm_difference, ]
  if (nrow(assays) == 0L)
    stop("No assays were found.", call. = FALSE)
  # Combine match matrices and calculate match for the entire assay
  assays <- combine_match_matrices(assays)
  assays <- dplyr::arrange(assays, begin)
  assays <- tibble::new_tibble(assays, class = "rprimer_assay")
  assays
}

#' Combine match matrices
#'
#' @param x A tibble with assays.
#'
#' @return A tibble with assays, with information on perfect matches   
#' for forward and reverse primers.
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
    match
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
  x
}

#' Add probes to match matrices
#'
#' @param x A tibble with assays (with probes).
#'
#' @return A tibble with assays, with information on perfect matches    
#' for forward and reverse primers and probes.
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
  x
}
