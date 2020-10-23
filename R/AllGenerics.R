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
#'
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
})

#' Plot an RprimerOligo-object (method)
#'
#' @param x An \code{RprimerOligo} object.
#'
#' @return A plot.
#'
#' @describeIn plotData
#'
#' @export
setMethod(
    "plotData", "RprimerOligo", function(x) {
        x <- as.data.frame(x)
        if (any(x$type == "probe")) {
            .primerProbePlot(x)
        } else {
            .primerPlot(x)
        }
})

#' Plot an RprimerAssay-object (method)
#'
#' @param x An \code{RprimerAssay} object.
#'
#' @return A plot.
#'
#' @describeIn plotData
#'
#' @export
setMethod(
    "plotData", "RprimerAssay", function(x) {
        x <- as.data.frame(x)
        .assayPlot(x)
})


# Helpers ======================================================================

.primerPlot <- function(x) {
    start <- x$alignmentStart[[1]]
    end <- x$alignmentEnd[[1]]
    colors <- c("#7B95A9", "#E7D0D8")
    ggplot2::ggplot() +
        .themeRprimer(showXAxis = TRUE) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
        ) +
        ggplot2::xlim(start, end) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "Position", y = "") +
        ggplot2::geom_segment(
            color = "grey", lwd = 2,
            ggplot2::aes(x = start, xend = end, y = 0, yend = 0)
        ) +
        ggplot2::geom_rect(
            data = x[x$type == "primer" & !is.na(x$majority), ], ggplot2::aes(
                xmin = start, xmax = end, ymin = 0.05 , ymax = 0.35
            ), fill = "#7B95A9"
        ) +
        ggplot2::geom_rect(
            data = x[x$type == "primer" & !is.na(x$majorityRc), ],
            ggplot2::aes(
                xmin = start, xmax = end, ymin = 0.35 , ymax = 0.65
            ), fill = "#E7D0D8"
        ) +
        ggplot2::annotate(
            "label", x = start, y = seq(0.89, length.out = 2, by = 0.07),
            label = c(
                "Forward primer", "Reverse primer"
            ), size = 3, hjust = 0, fontface = 2,
            color = colors, fill = "white", label.size = NA
        )
}

.primerProbePlot <- function(x) {
    start <- x$alignmentStart[[1]]
    end <- x$alignmentEnd[[1]]
    colors <- c("#7B95A9", "#E7D0D8", "#D09F99", "#404038")
    ggplot2::ggplot() +
        .themeRprimer(showXAxis = TRUE) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
        ) +
        ggplot2::xlim(start, end) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "Position", y = "") +
        ggplot2::geom_segment(
            color = "grey", lwd = 2,
            ggplot2::aes(x = start, xend = end, y = 0, yend = 0)
        ) +
        ggplot2::geom_rect(
            data = x[x$type == "primer" & !is.na(x$majority), ], ggplot2::aes(
            xmin = start, xmax = end, ymin = 0.05 , ymax = 0.20
            ), fill = "#7B95A9"
        ) +
        ggplot2::geom_rect(
            data = x[x$type == "primer" & !is.na(x$majorityRc), ],
            ggplot2::aes(
                xmin = start, xmax = end, ymin = 0.20 , ymax = 0.35
            ), fill = "#E7D0D8"
        ) +
        ggplot2::geom_rect(
            data = x[x$type == "probe" & !is.na(x$majority), ], ggplot2::aes(
                xmin = start, xmax = end, ymin = 0.35 , ymax = 0.5
            ), fill = "#D09F99"
        ) +
        ggplot2::geom_rect(
            data = x[x$type == "probe" & !is.na(x$majorityRc), ], ggplot2::aes(
                xmin = start, xmax = end, ymin = 0.5 , ymax = 0.65
            ), fill = "#404038"
        ) +
        ggplot2::annotate(
            "label", x = start, y = seq(0.75, length.out = 4, by = 0.07),
            label = c(
                "Forward primer", "Reverse primer", "Probe (+)", "Probe (-)"
            ), size = 3, hjust = 0, fontface = 2,
            color = colors, fill = "white", label.size = NA
        )
}

.assayPlot <- function(x) {
    start <- x$alignmentStart[[1]]
    end <- x$alignmentEnd[[1]]
    row <- seq_len(nrow(x))
    x <- tibble::add_column(x, row)
    colors <- c("#7B95A9", "#E7D0D8")
    ggplot2::ggplot() +
        .themeRprimer(showXAxis = TRUE) +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
        ) +
        ggplot2::xlim(start, end) +
        ggplot2::ylim(0, nrow(x)*1.05) +
        ggplot2::labs(x = "Position", y = "") +
        ggplot2::geom_segment(
            color = "grey", lwd = 2,
            ggplot2::aes(x = start, xend = end, y = 0, yend = 0)
        ) +
        ggplot2::geom_rect(
            data = x, ggplot2::aes(
                xmin = start, xmax = end, ymin = row, ymax = row + 0.3
            ), fill = "#7B95A9"
        ) +
        ggplot2::annotate(
            "label", x = start, y = nrow(x),
            label = paste(
                "Number of candiates:", nrow(x)
            ), size = 3, hjust = 0, fontface = 2,
            color = "#7B95A9", fill = "white", label.size = NA
        )
}
