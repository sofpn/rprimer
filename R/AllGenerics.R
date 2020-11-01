#' Plot an Rprimer-object (generic)
#'
#' @param x
#' An \code{RprimerProfile}, \code{RprimerOligo} or \code{RprimerAssay} object.
#'
#' @param ... Additional arguments that should be passed to the plot.
#'
#' @examples
#' data("exampleRprimerProfile")
#' plotData(exampleRprimerProfile)
#' plotData(exampleRprimerProfile, shadeFrom = 500, shadeTo = 1000)
#' data("exampleRprimerOligo")
#' plotData(exampleRprimerOligo)
#' data("exampleRprimerAssay")
#' plotData(exampleRprimerAssay)
#' @export
setGeneric("plotData", function(x, ...) standardGeneric("plotData"))
