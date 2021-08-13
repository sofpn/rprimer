#' rprimer Shiny application
#'
#' \code{runRprimerApp()} starts a Shiny application where
#' the workflow of the rprimer package can be run through a
#' graphical user interface.
#'
#' @return
#' Opens the Shiny application.
#'
#' @export
#'
#' @import shiny
#'
#' @examples
#' ## Only run this in interactive R sessions:
#' if (interactive()) {
#'     runRprimerApp()
#' }
runRprimerApp <- function() {
    rprimerApp()
}
