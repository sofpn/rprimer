#' Run the rprimer package as a shiny application
#'
#' @export
#'
#' @examples
#' runRprimerApp()
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
