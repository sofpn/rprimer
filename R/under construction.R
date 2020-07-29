




check_match_rprimer_assay <- function(x, y) {
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
  x <- dplyr::bind_cols(x, match_percentage)
  match_matrix <- tibble::tibble(match_matrix)
  x <- dplyr::bind_cols(x, match_matrix)
  x <- tibble::new_tibble(x, nrow = nrow(x), class = "rprimer_assay")
  return(x)
}

