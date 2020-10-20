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
#' @export
setGeneric("plotData", function(x, ...) standardGeneric("plotData"))

#' Plot an RprimerProfile-object (method)
#'
#' @param x An \code{RprimerProfile} object.
#'
#' @param shadeFrom
#' Optional. If a particular area should be shaded, where the shade
#' should start (an integer).
#'
#' @param shadeTo
#' Optional. If a particular area should be shaded, where the shade
#' should end (an integer).
#'
#' @return A plot.
#'
#' @describeIn plotData
#'
#' @importFrom patchwork wrap_plots
#'
#' @export
setMethod(
    "plotData", "RprimerProfile", function(x, shadeFrom = NULL, shadeTo = NULL) {
        if (is.null(shadeFrom) && is.null(shadeTo)) {
            shadeFrom <- -Inf
            shadeTo <- -Inf
        }
        if (is.null(shadeFrom) && !is.null(shadeTo)) {
            shadeFrom <- 0
        }
        if (!is.null(shadeFrom) && is.null(shadeTo)) {
            shadeTo <- nrow(x)
        }
        if (!is.null(shadeFrom) && !is.numeric(shadeFrom)) {
            stop("'shadeFrom' must be a number.", call. = FALSE)
        }
        if (!is.null(shadeTo) && !is.numeric(shadeTo)) {
            stop("'shadeTo' must be a number.", call. = FALSE)
        }
        if (ncol(x) != 11) {
            stop("Invalid number of columns for 'x'.", call. = FALSE)
        }
        x <- as.data.frame(x)
        patchwork::wrap_plots(
            list(
                .identityPlot(x, shadeFrom = shadeFrom, shadeTo = shadeTo),
                .entropyPlot(x, shadeFrom = shadeFrom, shadeTo = shadeTo),
                .gcPlot(x, shadeFrom = shadeFrom, shadeTo = shadeTo),
                .gapPlot(x, shadeFrom = shadeFrom, shadeTo = shadeTo)
            ),
            ncol = 1
        )
    }
)
