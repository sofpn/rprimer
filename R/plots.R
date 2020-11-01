# assay feature plot # also histograms/barchart for oligo features
# color blind palette

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
    }
)

#' Plot an RprimerOligo-object (method)
#'
#' @param x An \code{RprimerOligo} object.
#'
#' @return A plot.
#'
#' @describeIn plotData
#'
#' @export
setMethod("plotData", "RprimerOligo", function(x) {
        x <- as.data.frame(x)
        patchwork::wrap_plots(
            list(.oligoPlot(x), .oligoFeaturePlot(x)), ncol = 1
        )
    }
)

#' Plot an RprimerAssay-object (method)
#'
#' @param x An \code{RprimerAssay} object.
#'
#' @return A plot.
#'
#' @describeIn plotData
#'
#' @export
setMethod("plotData", "RprimerAssay", function(x) {
        x <- as.data.frame(x)
        patchwork::wrap_plots(
            list(.assayPlot(x), .assayFeaturePlot(x)), ncol = 1
        )
    }
)

#' Plot nucleotide distribution
#'
#' @param x An \code{RprimerProfile} object.
#'
#' @param rc
#' \code{TRUE} or {FALSE}. If the plotted sequence should be displayed
#' as reverse complement or not. Defaults to \code{FALSE}.
#'
#' @return A barplot with nucleotide distribution.
#'
#' @examples
#' data("exampleRprimerProfile")
#' plotNucleotides(exampleRprimerProfile[1:10, ])
#' ## Plot the first 10 bases
#' plotNucleotides(exampleRprimerProfile[1:10, ], rc = TRUE)
#' ## As reverse complement
#' @export
plotNucleotides <- function(x, rc = FALSE) {
    if (!methods::is(x, "RprimerProfile")) {
        stop("'x' must be an RprimerProfile object.", call. = FALSE)
    }
    if (!(is.logical(rc))) {
        stop("'rc' must be set to 'TRUE' or 'FALSE'.", call. = FALSE)
    }
    x <- as.data.frame(x)
    x <- dplyr::select(x, c("a", "c", "g", "t", "other"))
    x <- t(as.matrix(x))
    rownames(x)[-5] <- toupper(rownames(x))[-5]
    if (rc) {
        rownames(x) <- unname(rprimerGlobals$complementLookup[rownames(x)])
        x <- x[order(rownames(x)), ]
        if ("other" %in% rownames(x)) {
            swap <- c(which(rownames(x) == "other"), nrow(x))
            x[swap, ] <- x[rev(swap), ]
            rownames(x)[swap] <- rownames(x)[rev(swap)]
        }
    }
    Base <- Position <- Frequency <- NULL
    x <- reshape2::melt(x)
    names(x) <- c("Base", "Position", "Frequency")
    x$Base <- factor(x$Base)
    x$Position <- factor(x$Position)
    if (rc == TRUE) {
        x$Position <- factor(x$Position, levels = rev(levels(x$Position)))
    }
    cols <- c("#7B95A9", "#E7D0D8", "#D09F99", "#404038", "gray")
    ggplot2::ggplot(
        data = x, ggplot2::aes(x = Position, y = Frequency, fill = Base)
    ) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_manual(values = cols) +
        .themeRprimer(showLegend = TRUE)
}

# Helpers ======================================================================

.oligoPlot <- function(x) {
    if (any(x$type == "probe")) {
        .primerProbePlot(x)
    } else {
        .primerPlot(x)
    }
}

.primerPlot <- function(x) {
    start <- x$alignmentStart[[1]]
    end <- x$alignmentEnd[[1]]
    colors <- c("#7B95A9", "#E7D0D8")
    nPrimer <- nrow(x[x$type == "primer" & !is.na(x$majority), ])
    nPrimerRc <- nrow(x[x$type == "primer" & !is.na(x$majorityRc), ])
    ggplot2::ggplot() +
        ggplot2::xlim(start, end) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "Position", y = "") +
        ggplot2::geom_segment(
            color = "grey", lwd = 2,
            ggplot2::aes(x = start, xend = end, y = 0, yend = 0)
        ) +
        ggplot2::geom_rect(
            data = x[x$type == "primer" & !is.na(x$majority), ], ggplot2::aes(
                xmin = start, xmax = end, ymin = 0.05, ymax = 0.35
            ), fill = "#7B95A9"
        ) +
        ggplot2::geom_rect(
            data = x[x$type == "primer" & !is.na(x$majorityRc), ],
            ggplot2::aes(
                xmin = start, xmax = end, ymin = 0.35, ymax = 0.65
            ), fill = "#E7D0D8"
        ) +
        ggplot2::annotate(
            "label",
            x = start, y = seq(0.89, length.out = 2, by = 0.07),
            label = c(
                paste("Forward primer n =", nPrimer),
                paste("Reverse primer n =", nPrimerRc)
            ), size = 3, hjust = 0, fontface = 2,
            color = colors, fill = "white", label.size = NA
        ) +
        .themeRprimer(showYAxis = FALSE)
}

.primerProbePlot <- function(x) {
    start <- x$alignmentStart[[1]]
    end <- x$alignmentEnd[[1]]
    colors <- c("#7B95A9", "#E7D0D8", "#D09F99", "#404038")
    nPrimer <- nrow(x[x$type == "primer" & !is.na(x$majority), ])
    nPrimerRc <- nrow(x[x$type == "primer" & !is.na(x$majorityRc), ])
    nProbe <- nrow(x[x$type == "probe" & !is.na(x$majority), ])
    nProbeRc <- nrow(x[x$type == "probe" & !is.na(x$majorityRc), ])
    ggplot2::ggplot() +
        ggplot2::xlim(start, end) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "Position", y = "") +
        ggplot2::geom_segment(
            color = "grey", lwd = 2,
            ggplot2::aes(x = start, xend = end, y = 0, yend = 0)
        ) +
        ggplot2::geom_rect(
            data = x[x$type == "primer" & !is.na(x$majority), ], ggplot2::aes(
                xmin = start, xmax = end, ymin = 0.05, ymax = 0.20
            ), fill = "#7B95A9"
        ) +
        ggplot2::geom_rect(
            data = x[x$type == "primer" & !is.na(x$majorityRc), ],
            ggplot2::aes(
                xmin = start, xmax = end, ymin = 0.20, ymax = 0.35
            ), fill = "#E7D0D8"
        ) +
        ggplot2::geom_rect(
            data = x[x$type == "probe" & !is.na(x$majority), ], ggplot2::aes(
                xmin = start, xmax = end, ymin = 0.35, ymax = 0.5
            ), fill = "#D09F99"
        ) +
        ggplot2::geom_rect(
            data = x[x$type == "probe" & !is.na(x$majorityRc), ], ggplot2::aes(
                xmin = start, xmax = end, ymin = 0.5, ymax = 0.65
            ), fill = "#404038"
        ) +
        ggplot2::annotate(
            "label",
            x = start, y = seq(0.75, length.out = 4, by = 0.07),
            label = c(
                paste("Forward primer n =", nPrimer),
                paste("Reverse primer n =", nPrimerRc),
                paste("Probe (+) n =", nProbe),
                paste("Probe (-) n =", nProbeRc)
            ), size = 3, hjust = 0, fontface = 2,
            color = colors, fill = "white", label.size = NA
        ) +
        .themeRprimer(showYAxis = FALSE)
}

.violinPlot <- function(data, y, title = "", color = "#7B95A9") {
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

.gcTmIdentityPlot <- function(x, color = "grey60", type = "Primers") {
    patchwork::wrap_plots(
        list(
            .violinPlot(
                x, x$gcMajority, paste(type, "\nGC-content"), color = color
            ),
            .violinPlot(x, x$tmMajority, "\nTm", color = color),
            .violinPlot(x, x$identity, "\nIdentity", color = color)
        ),
        ncol = 3
    )
}

.oligoFeaturePlot <- function(x) {
    if (any(x$type == "probe")) {
        patchwork::wrap_plots(list(
            .gcTmIdentityPlot(x[x$type == "probe", ], type = "Probes"),
            .gcTmIdentityPlot(x[x$type == "primer", ], type = "Primers")
        ), ncol = 1)
    } else {
        patchwork::wrap_plots(list(
            .gcTmIdentityPlot(x[x$type == "primer", ], type = "Primers")
        ), ncol = 1)
    }
}

.assayPlot <- function(x) {
    start <- x$alignmentStart[[1]]
    end <- x$alignmentEnd[[1]]
    row <- seq_len(nrow(x))
    x <- tibble::add_column(x, row)
    colors <- c("#7B95A9", "#E7D0D8")
    ggplot2::ggplot() +
        ggplot2::xlim(start, end) +
        ggplot2::ylim(0, nrow(x) * 1.05) +
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
                "Assays, n =", nrow(x)
            ), size = 3, hjust = 0, fontface = 2,
            color = "#7B95A9", fill = "white", label.size = NA
        ) +
        .themeRprimer(showYAxis = FALSE)
}

.assayFeaturePlot <- function(x, color = "grey30") {
    patchwork::wrap_plots(
        list(
            .violinPlot(x, x$totalDegeneracy, "Total degeneracy", color = color),
            .violinPlot(x, x$meanIdentity, "Mean identity", color = color)
        ),
        ncol = 2
    )
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
            data = averages,  color = "#404038",
            ggplot2::aes(x = position, y = average)
        ) +
        ggplot2::ylim(0, 1) +
        ggplot2::ylab("Identity") +
        ggplot2::xlab("") +
        ggplot2::scale_color_gradient(low = "#a9bac785", high = "#a9bac785") +
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
            data = averages,  color = "#404038",
            ggplot2::aes(x = position, y = average)
        ) +
        ggplot2::ylab("Entropy") +
        ggplot2::xlab("") +
        ggplot2::scale_color_gradient(low = "#a9bac785", high = "#a9bac785") +
        .themeRprimer(showXAxis = FALSE)
}

.gcPlot <- function(x, shadeFrom = NULL, shadeTo = NULL) {
    position <- average <- NULL
    averages <- .gcRunningAverage(x$majority)
    xadj <- unique(x$position - seq_along(x$position))
    averages$position <- averages$position + xadj
    ggplot2::ggplot(data = x, ggplot2::aes(x = position)) +
        ggplot2::geom_segment(
            color = "#a9bac7",
            ggplot2::aes(
                x = min(position), xend = max(position), y = 0.5, yend = 0.5
            )
        ) +
        .shadeArea(shadeFrom = shadeFrom, shadeTo = shadeTo) +
        ggplot2::geom_line(
            data = averages, color = "#404038",
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
        data = x, ggplot2::aes(x = position, y = gaps)) +
        .shadeArea(shadeFrom = shadeFrom, shadeTo = shadeTo) +
        ggplot2::geom_point(size = 0.3, ggplot2::aes(color = gaps)) +
        ggplot2::ylim(0, 1) +
        ggplot2::xlab("Position") +
        ggplot2::ylab("Gaps") +
        ggplot2::scale_color_gradient(low = "#a9bac785", high = "#a9bac785") +
        .themeRprimer()
}

.shadeArea <- function(shadeFrom = NULL, shadeTo = NULL) {
    ggplot2::annotate(
        "rect",
        xmin = shadeFrom, xmax = shadeTo, ymin = -Inf, ymax = Inf,
        color = "grey80", alpha = 0.2
    )
}

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
