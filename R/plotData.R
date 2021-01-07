# else if
# restore runavg size...

#' Plot an Rprimer-object (generic)
#'
#' @param x
#' An \code{RprimerProfile}, \code{RprimerOligo} or \code{RprimerAssay} object.
#'
#' @param ... Additional arguments that should be passed to the plot.
#'
#' @return A plot.
#'
#' @examples
#' data("exampleRprimerProfile")
#' plotData(exampleRprimerProfile)
#' plotData(exampleRprimerProfile, shadeFrom = 500, shadeTo = 1000)
#' data("exampleRprimerOligo")
#' plotData(exampleRprimerOligo)
#' ## data("exampleRprimerAssay") ################ ADD LATER
#' ## plotData(exampleRprimerAssay)
#' @export
setGeneric("plotData", function(x, ...) standardGeneric("plotData"))

#' Plot an RprimerProfile-object (method)
#'
#' @param x An \code{RprimerProfile} object.
#'
#' @param shadeFrom
#' Optional. If a particular area should be shaded, where the shade
#' should begin (an integer).
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
setMethod("plotData", "RprimerProfile", function(x,
                                                 shadeFrom = NULL,
                                                 shadeTo = NULL) {
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
#' @importFrom patchwork wrap_plots
#'
#' @export
setMethod("plotData", "RprimerOligo", function(x) {
    x <- as.data.frame(x)
    patchwork::wrap_plots(
        list(.oligoPlot(x), .oligoFeaturePlot(x)),
        ncol = 1
    )
})

#' Plot an RprimerAssay-object (method)
#'
#' @param x An \code{RprimerAssay} object.
#'
#' @return A plot.
#'
#' @describeIn plotData
#'
#' @importFrom patchwork wrap_plots
#'
#' @export
setMethod("plotData", "RprimerAssay", function(x) {
    x <- as.data.frame(x)
    patchwork::wrap_plots(
        list(.assayPlot(x), .assayFeaturePlot(x)),
        ncol = 1
    )
})


# Helpers ======================================================================


.addEmptyRow <- function(x) {
    x[1, ] <- NA
    x$start <- 0
    x$end <- 0
    x
}

.splitOligoDf <- function(x) {
    all <- list()
    all$fwd <- x[x$type == "primer" & x$fwd, , drop = FALSE]
    all$rev <- x[x$type == "primer" & x$rev, , drop = FALSE]
    if (any(x$type == "probe")) {
        all$prPos <- x[x$type == "probe" & x$fwd, , drop = FALSE]
        all$prNeg <- x[x$type == "probe" & x$rev, , drop = FALSE]
    }
    nRows <- vapply(all, nrow, integer(1), USE.NAMES = FALSE)
    all[nRows == 0] <- lapply(all[nRows == 0], .addEmptyRow)
    all
}

.oligoPlot <- function(x) {
    if (all(x$type == "primer")) {
        .primerPlot(x)
    } else if (all(x$type == "probe")) {
        .probePlot(x)
    } else {
        .primerProbePlot(x)
    }
}

.primerPlot <- function(x) {
    start <- end <- NULL
    roiStart <- x$roiStart[[1]]
    roiEnd <- x$roiEnd[[1]]
    colors <- c(fwd = "#A4A7B6", rev = "#82879B")
    oligos <- .splitOligoDf(x)
    ggplot2::ggplot() +
        ggplot2::xlim(roiStart, roiEnd) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "Position", y = "") +
        ggplot2::geom_segment(
            color = "grey", lwd = 2, ggplot2::aes(
                x = roiStart, xend = roiEnd, y = 0, yend = 0
            )
        ) +
        ggplot2::geom_rect(data = oligos$fwd, ggplot2::aes(
            xmin = start, xmax = end, ymin = 0.35, ymax = 0.65
        ), fill = colors["fwd"]) +
        ggplot2::geom_rect(data = oligos$rev, ggplot2::aes(
            xmin = start, xmax = end, ymin = 0.05, ymax = 0.35
        ), fill = colors["rev"]) +
        ggplot2::annotate(
            "label",
            x = roiStart,
            y = seq(0.89, length.out = 2, by = 0.07), label = c(
                paste(
                    "Reverse primer n =",
                    nrow(oligos$rev[!is.na(oligos$rev$length), ])
                ),
                paste(
                    "Forward primer n =",
                    nrow(oligos$fwd[!is.na(oligos$fwd$length), ])
                )
            ), size = 3, hjust = 0, fontface = 2,
            color = rev(colors), fill = "white", label.size = NA
        ) +
        .themeRprimer(showYAxis = FALSE)
}

.probePlot <- function(x) {
    start <- end <- NULL
    roiStart <- x$roiStart[[1]]
    roiEnd <- x$roiEnd[[1]]
    colors <- c(prPos = "#64697D", prNeg = "#525666")
    oligos <- .splitOligoDf(x)
    ggplot2::ggplot() +
        ggplot2::xlim(roiStart, roiEnd) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "Position", y = "") +
        ggplot2::geom_segment(
            color = "grey", lwd = 2, ggplot2::aes(
                x = roiStart, xend = roiEnd, y = 0, yend = 0
            )
        ) +
        ggplot2::geom_rect(data = oligos$prPos, ggplot2::aes(
            xmin = start, xmax = end, ymin = 0.35, ymax = 0.65
        ), fill = colors["prPos"]) +
        ggplot2::geom_rect(data = oligos$prNeg, ggplot2::aes(
            xmin = start, xmax = end, ymin = 0.05, ymax = 0.35
        ), fill = colors["prNeg"]) +
        ggplot2::annotate(
            "label",
            x = roiStart,
            y = seq(0.89, length.out = 2, by = 0.07), label = c(
                paste(
                    "Probe (-) n =",
                    nrow(oligos$prNeg[!is.na(oligos$prNeg$length), ])
                ),
                paste(
                    "Probe (+) n =",
                    nrow(oligos$prPos[!is.na(oligos$prPos$length), ])
                )
            ), size = 3, hjust = 0, fontface = 2,
            color = rev(colors), fill = "white", label.size = NA
        ) +
        .themeRprimer(showYAxis = FALSE)
}

.primerProbePlot <- function(x) {
    start <- end <- NULL
    roiStart <- x$roiStart[[1]]
    roiEnd <- x$roiEnd[[1]]
    oligos <- .splitOligoDf(x)
    colors <- c(
        fwd = "#A4A7B6", rev = "#82879B", prPos = "#64697D", prNeg = "#525666"
    )
    ggplot2::ggplot() +
        ggplot2::xlim(roiStart, roiEnd) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "Position", y = "") +
        ggplot2::geom_segment(
            color = "grey", lwd = 2, ggplot2::aes(
                x = roiStart, xend = roiEnd, y = 0, yend = 0
            )
        ) +
        ggplot2::geom_rect(data = oligos$fwd, ggplot2::aes(
            xmin = start, xmax = end, ymin = 0.5, ymax = 0.65
        ), fill = colors["fwd"]) +
        ggplot2::geom_rect(data = oligos$rev, ggplot2::aes(
            xmin = start, xmax = end, ymin = 0.35, ymax = 0.5
        ), fill = colors["rev"]) +
        ggplot2::geom_rect(data = oligos$prPos, ggplot2::aes(
            xmin = start, xmax = end, ymin = 0.20, ymax = 0.35
        ), fill = colors["prPos"]) +
        ggplot2::geom_rect(data = oligos$prNeg, ggplot2::aes(
            xmin = start, xmax = end, ymin = 0.05, ymax = 0.20
        ), fill = colors["prNeg"]) +
        ggplot2::annotate(
            "label",
            x = roiStart, y = seq(0.75, length.out = 4, by = 0.07),
            label = c(
                paste(
                    "Probe (-) n =",
                    nrow(oligos$prNeg[!is.na(oligos$prNeg$length), ])
                ),
                paste(
                    "Probe (+) n =",
                    nrow(oligos$prPos[!is.na(oligos$prPos$length), ])
                ),
                paste(
                    "Reverse primer n =",
                    nrow(oligos$rev[!is.na(oligos$rev$length), ])
                ),
                paste(
                    "Forward primer n =",
                    nrow(oligos$fwd[!is.na(oligos$fwd$length), ])
                )
            ), size = 3, hjust = 0, fontface = 2,
            color = rev(colors), fill = "white", label.size = NA
        ) +
        .themeRprimer(showYAxis = FALSE)
}

.violinPlot <- function(data, y, title = "", color = "grey30") {
    ggplot2::ggplot() +
        ggplot2::geom_violin(
            data = data, ggplot2::aes(x = 1, y = y),
            fill = color, color = NA, alpha = 0.4
        ) +
        ggplot2::geom_point(
            data = data, ggplot2::aes(x = 1, y = y), alpha = 0.05, color = color
        ) +
        ggplot2::geom_boxplot(
            data = data, ggplot2::aes(x = 1, y = y), width = 0.1,
            color = color, fill = color, alpha = 0.2
        ) +
        ggplot2::ylab("") +
        ggplot2::labs(title = title) +
        .themeRprimer(showXAxis = FALSE)
}

.boxPlot <- function(data, y, title = "", color = "grey30") {
    ggplot2::ggplot() +
        ggplot2::geom_point(
            data = data, ggplot2::aes(x = 1, y = y), alpha = 0.05, color = color
        ) +
        ggplot2::geom_boxplot(
            data = data, ggplot2::aes(x = 1, y = y), width = 0.1,
            color = color, fill = color, alpha = 0.2
        ) +
        ggplot2::ylab("") +
        ggplot2::labs(title = title) +
        .themeRprimer(showXAxis = FALSE)
}

.dotPlot <- function(data, y, title = "", color = "grey30") {
    ggplot2::ggplot() +
        ggplot2::geom_point(
            data = data, ggplot2::aes(x = 1, y = y), alpha = 0.5, color = color
        ) +
        ggplot2::ylab("") +
        ggplot2::labs(title = title) +
        .themeRprimer(showXAxis = FALSE)
}

.barPlot <- function(data, y, title = "", color = "grey30") {
    ggplot2::ggplot(data, ggplot2::aes(factor(y))) +
        ggplot2::geom_bar(
            fill = color, color = color, alpha = 0.4
        ) +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::labs(title = title) +
        .themeRprimer(showXAxis = TRUE)
}

.gcTmIdentityPlot <- function(x, color = "grey30", type = "Primers") {
    if (nrow(x) >= 10) {
        patchwork::wrap_plots(
            list(
                .violinPlot(
                    x, x$gcContentMean,
                    paste0("\n", type, "\n\nGC-content"),
                    color = color
                ),
                .violinPlot(x, x$tmMean, "\n\n\nTm", color = color),
                .violinPlot(x, x$identity, "\n\n\nIdentity", color = color),
                .barPlot(x, x$length, "\n\n\nLength", color = color),
                .barPlot(x, x$degeneracy, "\n\n\nDegeneracy", color = color)
            ),
            ncol = 5
        )
    } else {
        patchwork::wrap_plots(
            list(
                .dotPlot(
                    x, x$gcContentMean,
                    paste0("\n", type, "\n\nGC-content"),
                    color = color
                ),
                .dotPlot(x, x$tmMean, "\n\n\nTm", color = color),
                .dotPlot(x, x$identity, "\n\n\nIdentity", color = color),
                .barPlot(x, x$length, "\n\n\nLength", color = color),
                .barPlot(x, x$degeneracy, "\n\n\nDegeneracy", color = color)
            ),
            ncol = 5
        )
    }
}

.oligoFeaturePlot <- function(x) {
    if (any(x$type == "probe")) {
        patchwork::wrap_plots(list(
            .gcTmIdentityPlot(
                x[x$type == "primer", ],
                color = "grey30", type = "Primers"
            ),
            .gcTmIdentityPlot(
                x[x$type == "probe", ],
                color = "grey30", type = "Probes"
            )
        ), ncol = 1)
    } else {
        patchwork::wrap_plots(list(
            .gcTmIdentityPlot(
                x[x$type == "primer", ],
                color = "grey30",
                type = "Primers"
            )
        ), ncol = 1)
    }
}

.assayPlot <- function(x) {
    start <- x$roiStart[[1]]
    end <- x$roiEnd[[1]]
    ggplot2::ggplot() +
        ggplot2::xlim(start, end) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "Position", y = "") +
        ggplot2::geom_segment(
            color = "grey", lwd = 2,
            ggplot2::aes(x = start, xend = end, y = 0, yend = 0)
        ) +
        ggplot2::geom_rect(
            data = x, ggplot2::aes(
                xmin = start, xmax = end, ymin = 0.05, ymax = 0.65
            ), fill = "#64697D"
        ) +
        ggplot2::annotate(
            "label",
            x = start, y = 0.8,
            label = paste(
                "Assays n =", nrow(x)
            ), size = 3, hjust = 0, fontface = 2,
            color = "#64697D", fill = "white", label.size = NA
        ) +
        .themeRprimer(showYAxis = FALSE)
}

.assayFeaturePlot <- function(x, color = "grey30") {
    if (nrow(x) >= 10) {
        patchwork::wrap_plots(
            list(
                .violinPlot(
                    x, x$tmDifferencePrimer, "\n\nTm differende, primers",
                    color = color
                ),
                .boxPlot(
                    x, x$ampliconLength, "\n\nAmplicon length",
                    color = color
                ),
                .barPlot(
                    x, x$totalDegeneracy, "\n\nTotal degeneracy",
                    color = color
                )
            ),
            ncol = 3
        )
    } else {
        patchwork::wrap_plots(
            list(
                .dotPlot(
                    x, x$tmDifferencePrimer, "\n\nTm difference, primers",
                    color = color
                ),
                .barPlot(
                    x, x$ampliconLength, "\n\nAmplicon length",
                    color = color
                ),
                .barPlot(
                    x, x$totalDegeneracy, "\n\nTotal degeneracy",
                    color = color
                )
            ),
            ncol = 3
        )
    }
}

# Helpers for plotting an RprimerProfile =======================================

#' Calculate running average
#'
#' @param x A numeric vector.
#'
#' @param size
#' The number of observations in each average.
#' If \code{NULL}, the size will be set to the nearest positive, nonzero
#' integer to \code{length(x)/100}.
#'
#' @return A data frame with position and running average of \code{x}.
#'
#' @keywords internal
#'
#' @noRd
.runningAverage <- function(x, size = NULL) {
    if (is.null(size)) {
        size <- round(length(x) / 100)
        if (size == 0) {
            size <- 1
        }
    }
    sums <- c(0, cumsum(x))
    from <- seq_len(length(sums) - size)
    to <- seq(size + 1, length(sums))
    average <- (sums[to] - sums[from]) / size
    position <- seq(size, length(x))
    if (size > 1) {
        midpoint <- size / 2
        position <- position - midpoint
    }
    data.frame(position, average)
}

.identityPlot <- function(x, shadeFrom = NULL, shadeTo = NULL) {
    position <- identity <- average <- NULL
    averages <- .runningAverage(x$identity)
    xadj <- unique(x$position - seq_along(x$position))
    averages$position <- averages$position + xadj
    ggplot2::ggplot(
        data = x, ggplot2::aes(x = position, y = identity)
    ) +
        .shadeArea(shadeFrom = shadeFrom, shadeTo = shadeTo) +
        ggplot2::geom_point(size = 0.3, ggplot2::aes(colour = identity)) +
        ggplot2::geom_line(
            data = averages,  color = "#1B1C22",
            ggplot2::aes(x = position, y = average)
        ) +
        ggplot2::ylim(0, 1) +
        ggplot2::ylab("Identity") +
        ggplot2::xlab("") +
        ggplot2::scale_color_gradient(low = "#BABFD1", high = "#BABFD1") +
        .themeRprimer(showXAxis = FALSE)
}

.entropyPlot <- function(x, shadeFrom = NULL, shadeTo = NULL) {
    position <- entropy <- average <- NULL
    averages <- .runningAverage(x$entropy)
    xadj <- unique(x$position - seq_along(x$position))
    averages$position <- averages$position + xadj
    ggplot2::ggplot(
        data = x, ggplot2::aes(x = position, y = entropy)
    ) +
        .shadeArea(shadeFrom = shadeFrom, shadeTo = shadeTo) +
        ggplot2::geom_point(size = 0.3, ggplot2::aes(colour = entropy)) +
        ggplot2::geom_line(
            data = averages,  color = "#1B1C22",
            ggplot2::aes(x = position, y = average)
        ) +
        ggplot2::ylab("Entropy") +
        ggplot2::xlab("") +
        ggplot2::scale_color_gradient(low = "#BABFD1", high = "#BABFD1") +
        .themeRprimer(showXAxis = FALSE)
}

.gcPlot <- function(x, shadeFrom = NULL, shadeTo = NULL) {
    position <- average <- NULL
    gc <- ifelse(x$majority == "C" | x$majority == "G", 1, 0)
    averages <- .runningAverage(gc)
    xadj <- unique(x$position - seq_along(x$position))
    averages$position <- averages$position + xadj
    ggplot2::ggplot(data = x, ggplot2::aes(x = position)) +
        ggplot2::geom_segment(
            color = "#BABFD1",
            ggplot2::aes(
                x = min(position), xend = max(position), y = 0.5, yend = 0.5
            )
        ) +
        .shadeArea(shadeFrom = shadeFrom, shadeTo = shadeTo) +
        ggplot2::geom_line(
            data = averages, color = "#1B1C22",
            ggplot2::aes(x = position, y = average)
        ) +
        ggplot2::xlab("") +
        ggplot2::ylab("GC-content") +
        ggplot2::ylim(0, 1) +
        .themeRprimer(showXAxis = FALSE)
}

.gapPlot <- function(x, shadeFrom = NULL, shadeTo = NULL) {
    position <- gaps <- NULL
    ggplot2::ggplot(
        data = x, ggplot2::aes(x = position, y = gaps)
    ) +
        .shadeArea(shadeFrom = shadeFrom, shadeTo = shadeTo) +
        ggplot2::geom_point(size = 0.3, ggplot2::aes(color = gaps)) +
        ggplot2::ylim(0, 1) +
        ggplot2::xlab("Position") +
        ggplot2::ylab("Gaps") +
        ggplot2::scale_color_gradient(low = "#BABFD1", high = "#BABFD1") +
        .themeRprimer()
}

.shadeArea <- function(shadeFrom = NULL, shadeTo = NULL) {
    ggplot2::annotate(
        "rect",
        xmin = shadeFrom, xmax = shadeTo, ymin = -Inf, ymax = Inf,
        color = "white", alpha = 0.2, fill = "grey30"
    )
}

# Theme ========================================================================

.themeRprimer <- function(showXAxis = TRUE,
                          showYAxis = TRUE,
                          showLegend = FALSE) {
    if (showXAxis && showYAxis) {
        .showXYAxes(showLegend = showLegend)
    } else if (showXAxis && !showYAxis) {
        .showXAxisHideYAxis(showLegend = showLegend)
    } else if (!showXAxis && showYAxis) {
        .showYAxisHideXAxis(showLegend = showLegend)
    } else {
        .hideXYAxes(showLegend = showLegend)
    }
}

.showXYAxes <- function(showLegend = TRUE) {
    ggplot2::theme_light() +
        ggplot2::theme(
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(color = "grey30", size = 9),
            legend.position = ifelse(showLegend, "right", "none"),
            axis.title.y = ggplot2::element_text(
                size = 9, color = "grey30",
                margin = ggplot2::margin(r = 10)
            ),
            axis.title.x = ggplot2::element_text(
                size = 9, color = "grey30",
                margin = ggplot2::margin(t = 10)
            ),
            plot.title = ggplot2::element_text(size = 9, color = "grey30"),
            plot.margin = ggplot2::unit(rep(0.1, 4), "cm")
        )
}

.showXAxisHideYAxis <- function(showLegend = TRUE) {
    ggplot2::theme_light() +
        ggplot2::theme(
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(color = "grey30", size = 9),
            legend.position = ifelse(showLegend, "right", "none"),
            axis.title.x = ggplot2::element_text(
                size = 9, color = "grey30",
                margin = ggplot2::margin(t = 10)
            ),
            axis.title.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(size = 9, color = "grey30"),
            plot.margin = ggplot2::unit(rep(0.1, 4), "cm")
        )
}

.showYAxisHideXAxis <- function(showLegend = TRUE) {
    ggplot2::theme_light() +
        ggplot2::theme(
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(color = "grey30", size = 9),
            legend.position = ifelse(showLegend, "right", "none"),
            axis.title.y = ggplot2::element_text(
                size = 9, color = "grey30",
                margin = ggplot2::margin(r = 10)
            ),
            axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(size = 9, color = "grey30"),
            plot.margin = ggplot2::unit(rep(0.1, 4), "cm")
        )
}

.hideXYAxes <- function(showLegend = TRUE) {
    ggplot2::theme_light() +
        ggplot2::theme(
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(color = "grey30", size = 9),
            legend.position = ifelse(showLegend, "right", "none"),
            axis.title.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(size = 9, color = "grey30"),
            plot.margin = ggplot2::unit(rep(0.1, 4), "cm")
        )
}
