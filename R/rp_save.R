#' Save an rprimer object to a file (generic)
#'
#' @param x An rprimer object (see methods for details).
#'
#' @param filename The name of the file. A character vector of length one.
#'
#' @details
#' The file format will depend on object.
#' Objects of class rprimer_alignment will be saved as .txt-files,
#' whereas objects of class rprimer_properties, rprimer_oligo
#' and rprimer_assay will be saved as .csv-files.
#'
#' @return A saved file.
#'
#' @examples
#' \dontrun{
#'   rp_save(example_rprimer_alignment, 'my_alignment')
#'   rp_save(example_rprimer_properties, 'my_sequence_properties')
#'   rp_save(example_rprimer_oligo, 'my_oligos')
#'   rp_save(example_rprimer_assay, 'my_assays')
#' }
#'
#' @export
rp_save <- function(x, filename) {
  if (!is.character(filename)) {
    stop(
      "'filename' must be a character vector.", call. = FALSE
    )
  }
  object_name <- as.character(substitute(x))
  if (!exists(object_name)) {
    stop(paste("object", object_name, "does not exist."), call. = FALSE)
  }
  UseMethod("rp_save")
}

#' @describeIn rp_save Save an rprimer_alignment object.
#'
#' @export
rp_save.rprimer_alignment <- function(x, filename) {
  seq_names <- purrr::map_chr(names(x), ~paste0(">", .x))
  my_file <- file(paste0(filename, ".txt"), open = "w")
  purrr::walk2(seq_names, x, ~writeLines(text = c(.x, .y), con = my_file))
  close(my_file)
}

#' @describeIn rp_save Save an rprimer_properties object.
#'
#' @export
rp_save.rprimer_properties <- function(x, filename) {
  utils::write.csv(
    x, file = paste0(filename, ".csv"), quote = FALSE, row.names = FALSE
  )
}

#' @describeIn rp_save Save an rprimer_oligo object.
#'
#' @export
rp_save.rprimer_oligo <- function(x, filename) {
  if (any(grepl("match_matrix", names(x)))) {
    x <- x[, names(x) != "match_matrix"]
  }
  utils::write.csv(
    x, file = paste0(filename, ".csv"), quote = FALSE, row.names = FALSE
  )
}

#' @describeIn rp_save Save an'rprimer_assay' object.
#'
#' @export
rp_save.rprimer_assay <- function(x, filename) {
  if (any(grepl("match_matrix", names(x)))) {
    x <- x[, names(x) != "match_matrix"]
  }
  utils::write.csv(
    x, file = paste0(filename, ".csv"), quote = FALSE, row.names = FALSE
  )
}
