#' rprimer Shiny application
#'
#' \code{runRprimerApp()} starts a Shiny application where
#' the workflow of the rprimer package can be run through a
#' graphical user interface.
#'
#' @return
#' Opens the Shiny application.
#'
#' @import shiny shinydashboard shinycssloaders
#'
#' @export
#'
#' @examples
#' ## Only run the application in interactive R sessions:
#' if (interactive()) {
#'     runRprimerApp()
#' }
runRprimerApp <- function() {
    appDir <- system.file("shiny-app", "app.R", package = "rprimer")
    if (appDir == "") {
        stop(
            "Could not find the app directory. Try re-installing `rprimer`.",
            call. = FALSE
        )
    }
    shiny::runApp(appDir, display.mode = "normal")
}
