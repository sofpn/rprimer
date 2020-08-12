#' Write assay report
#'
#' @param filename (a character vector of length one).
#'
#' @param assay_selection An object of class 'rprimer_assay', with one row.
#'
#' @param sequence_profile An object of class 'rprimer_profile'.
#'
#' @param sequence_properties An object of class 'rprimer_properties.
#'
#' @param comment
#' Optional. Comments that should be included in the report.
#' A character vector of length one. The default is \code{NULL}.
#'
#' @return A html-document with detailed assay information.
#'
#' @details Details:
#' rmarkdown and kableExtra are needed for this function to work.
#'
#' @export
write_report <- function(filename = "my_assay_report",
                         assay_selection,
                         sequence_profile,
                         sequence_properties,
                         comment = NULL) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop(
      "rmkardown is needed for this function to work. \n
      Please install it.", call. = FALSE
    )
  }
  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    stop(
      "kableExtra is needed for this function to work. \n
      Please install it.", call. = FALSE
    )
  }
  if (typeof(filename) != "character" || length(filename) != 1) {
    stop("'filename' must be a character vector.", call. = FALSE)
  }
  if (!inherits(assay_selection, "rprimer_assay")) {
    stop(
      "'assay_selection' must be an rprimer_assay object.", #### If several rows?
      call. = FALSE
    )
  }
  if (!is.rprimer_profile(sequence_profile)) {
    stop(
      "'sequence_profile' must be an rprimer_profile object.",
      call. = FALSE
    )
  }
  if (!is.rprimer_properties(sequence_properties)) {
    stop(
      "'sequence_properties' must be an rprimer_properties object.",
      call. = FALSE
    )
    if (is.null(comment)) comment <- ""
    if (!is.character(comment)) {
      stop("'comment' must be a character vector.")
    }
  }
  filename <- paste0(filename, ".html")
  rmarkdown::render(
    input = "assay_report.Rmd",
    output_file = filename,
    params = list(
      assay_selection = assay_selection,
      sequence_profile = sequence_profile,
      sequence_properties = sequence_properties,
      comment = comment
    )
  )
}
