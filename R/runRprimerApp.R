#' rprimer Shiny application
#'
#' \code{runRprimerApp()} starts a Shiny application where
#' the workflow of the rprimer package can be run through a
#' graphical user interface.
#'
#' @return
#' Opens the Shiny application.
#'
#' @import bslib shiny shinycssloaders shinyFeedback
#'
#' @export
#'
#' @examples
#' ## Only run this in interactive R sessions:
#' if (interactive()) {
#'     runRprimerApp()
#' }
runRprimerApp <- function() {
    dir <- system.file("shiny-app", "app.R", package = "rprimer")
    if (dir == "") {
        stop(
            "The app directory was not found. Try to install `rprimer` again.",
            call. = FALSE
        )
    }
    shiny::runApp(dir, display.mode = "normal")
}
