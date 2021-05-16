#' Plot an Rprimer object (generic)
#'
#' @param x
#' An \code{RprimerProfile}, \code{RprimerOligo} or \code{RprimerAssay} object.
#'
#' @param ...
#' Optional arguments for \code{RprimerProfile} objects.
#'
#' @param type
#' For \code{Rprimeroligo} objects:
#' Type of plot: \code{"overview"}, or
#' \code{"nucleotide"}, defaults to \code{"overview"}.
#'
#' @param highlight
#' For \code{Rprimeroligo} objects:
#' if a specific region within an overview plot should be highlighted.
#' A numeric vector indicating the start and end position,
#' e.g. \code{c(100, 1000)}, defaults to \code{NULL}
#' (i.e., no highlight).
#'
#' @param rc
#' For \code{Rprimeroligo} objects, and \code{type = "nucleotide"}:
#' If the plotted sequence should be displayed
#' as reverse complement or not.
#' \code{TRUE} or {FALSE}, defaults to \code{FALSE}.
#'
#' See examples below.
#'
#' @return A plot.
#'
#' @export
setGeneric("plotData", function(x, ...) standardGeneric("plotData"))

# Methods ======================================================================

#' Plot an RprimerProfile object (method)
#'
#' @describeIn plotData
#'
#' @importFrom patchwork wrap_plots
#'
#' @export
#'
#' @examples
#' #### Plot an RprimerProfile object
#'
#' data("exampleRprimerProfile")
#' prof <- exampleRprimerProfile
#'
#' ## Plot an overwiev
#' plotData(prof)
#'
#' ## Highlight a specific area
#' plotData(prof, highlight = c(500, 1000))
#'
#' ## Select a region of interest
#' roi <- prof[prof$position >= 500 & prof$position <= 550, ]
#'
#' ## Plot an overview of the roi
#' plotData(roi)
#'
#' ## Plot the nucleotide distribution of the roi
#' plotData(roi, type = "nucleotide")
#'
#' ## Plot the nucleotide distribution of the roi, as reverse complement
#' plotData(roi, type = "nucleotide", rc = TRUE)
setMethod("plotData", "RprimerProfile", function(x,
                                                 type = "overview",
                                                 highlight = NULL,
                                                 rc = FALSE) {
    if (nrow(x) == 0L) {
        stop("'x' does not contain any observations.", call. = FALSE)
    }
    if (type == "overview") {
        if (is.null(highlight)) {
            highlight <- -Inf
        }
        if (!is.numeric(highlight)) {
            stop("'highlight' must be numeric, e.g. c(1, 10).", call. = FALSE)
        }
        .plotOverview(x, highlight)
    }
    else if (type == "nucleotide") {
        if (!(is.logical(rc))) {
            stop("'rc' must be set to TRUE or FALSE.", call. = FALSE)
        }
        .plotNucleotides(x, rc)
    }
    else {
        stop("'type' must be either 'overview' or 'nucleotide'", call. = FALSE)
    }
})

#' Plot an RprimerOligo object (method)
#'
#' @describeIn plotData
#'
#' @aliases plotData
#'
#' @importFrom patchwork wrap_plots
#'
#' @export
#'
#' @examples
#' ### Plot an RprimerOligo object
#'
#' data("exampleRprimerOligo")
#' oligos <- exampleRprimerOligo
#'
#' ## Plot all oligos
#' plotData(oligos)
#'
#' ## Select a subset of the oligos, and plot
#' selected <- oligos[oligos$start >= 5000, ]
#' plotData(selected)
setMethod("plotData", "RprimerOligo", function(x) {
    if (nrow(x) == 0L) {
        stop("'x' does not contain any observations.", call. = FALSE)
    }
    x <- as.data.frame(x)
    patchwork::wrap_plots(
        list(.oligoPlot(x), .oligoFeaturePlot(x)),
        ncol = 1
    )
})

#' Plot an RprimerAssay object (method)
#'
#' @describeIn plotData
#'
#' @importFrom patchwork wrap_plots
#'
#' @export
#'
#' @examples
#' ### Plot an RprimerAssay object
#'
#' data("exampleRprimerAssay")
#' plotData(exampleRprimerAssay)
setMethod("plotData", "RprimerAssay", function(x) {
    if (nrow(x) == 0L) {
        stop("'x' does not contain any observations.", call. = FALSE)
    }
    x <- as.data.frame(x)
    patchwork::wrap_plots(
        list(.assayPlot(x), .assayFeaturePlot(x)),
        ncol = 1
    )
})

#' Plot an RprimerMatchOligo object (method)
#'
#' @describeIn plotData
#'
#' @export
#'
#' @examples
#' ### Plot an RprimerMatchOligo object
#'
#' data("exampleRprimerMatchOligo")
#' plotData(exampleRprimerMatchOligo)
setMethod("plotData", "RprimerMatchOligo", function(x) {
    if (nrow(x) == 0L) {
        stop("'x' does not contain any observations.", call. = FALSE)
    }
    x <- as.data.frame(x)
    .plotMatch(x)
})

#' Plot an RprimerMatchAssay object (method)
#'
#' @describeIn plotData
#'
#' @export
#'
#' @examples
#' ### Plot an RprimerMatchAssay object
#'
#' data("exampleRprimerMatchAssay")
#' plotData(exampleRprimerMatchAssay)
setMethod("plotData", "RprimerMatchAssay", function(x) {
    if (nrow(x) == 0L) {
        stop("'x' does not contain any observations.", call. = FALSE)
    }
    x <- as.data.frame(x)
    .plotMatchAssay(x)
})

# Helpers for plotting an RprimerOligo/Assay ===================================

.addEmptyRow <- function(x) {
    x[1, ] <- NA
    x$start <- 1
    x$end <- 1
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
    colors <- c(fwd = "#424B54", rev = "#93A8AC")
    oligos <- .splitOligoDf(x)
    ggplot2::ggplot() +
        ggplot2::xlim(roiStart, roiEnd) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "Position", y = "") +
        ggplot2::geom_segment(
            color = "grey60", lwd = 2, ggplot2::aes(
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
    colors <- c(prPos = "#9B6A6C", prNeg = "#E2B4BD")
    oligos <- .splitOligoDf(x)
    ggplot2::ggplot() +
        ggplot2::xlim(roiStart, roiEnd) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "Position", y = "") +
        ggplot2::geom_segment(
            color = "grey60", lwd = 2, ggplot2::aes(
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
        fwd = "#424B54", rev = "#93A8AC", prPos = "#9B6A6C", prNeg = "#E2B4BD"
    )
    ggplot2::ggplot() +
        ggplot2::xlim(roiStart, roiEnd) +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(x = "Position", y = "") +
        ggplot2::geom_segment(
            color = "grey60", lwd = 2, ggplot2::aes(
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

.violinPlot <- function(data, y, title = "", color = "grey20") {
    ggplot2::ggplot() +
        ggplot2::geom_violin(
            data = data, ggplot2::aes(x = 1, y = y),
            fill = color, alpha = 0.4
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

.boxPlot <- function(data, y, title = "", color = "grey20") {
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

.dotPlot <- function(data, y, title = "", color = "grey20") {
    ggplot2::ggplot() +
        ggplot2::geom_point(
            data = data, ggplot2::aes(x = 1, y = y), alpha = 0.5, color = color
        ) +
        ggplot2::ylab("") +
        ggplot2::labs(title = title) +
        .themeRprimer(showXAxis = FALSE)
}

.barPlot <- function(data, y, title = "", color = "grey20") {
    ggplot2::ggplot(data, ggplot2::aes(factor(y))) +
        ggplot2::geom_bar(
            fill = color, color = color, alpha = 0.4
        ) +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::labs(title = title) +
        .themeRprimer(showXAxis = TRUE)
}

.gcTmCoveragePlot <- function(x, color = "grey20", type = "Primers") {
    if (nrow(x) >= 10) {
        patchwork::wrap_plots(
            list(
                .violinPlot(
                    x, x$gcContentMean,
                    paste0("\n", type, "\n\nGC-content (mean)"),
                    color = color
                ),
                .violinPlot(x, x$tmMean, "\n\n\nTm (mean)", color = color),
                .violinPlot(
                    x, x$identity, "\n\n\nIdentity (mean)",
                    color = color
                ),
                .violinPlot(
                    x, x$coverage, "\n\n\nCoverage (mean)",
                    color = color
                ),
                .barPlot(x, x$length, "\n\n\nLength", color = color),
                .barPlot(x, x$degeneracy, "\n\n\nDegeneracy", color = color)
            ),
            ncol = 6
        )
    } else {
        patchwork::wrap_plots(
            list(
                .dotPlot(
                    x, x$gcContentMean,
                    paste0("\n", type, "\n\nGC-content (mean)"),
                    color = color
                ),
                .dotPlot(x, x$tmMean, "\n\n\nTm (mean)", color = color),
                .dotPlot(
                    x, x$identity, "\n\n\nIdentity (mean)",
                    color = color
                ),
                .dotPlot(
                    x, x$coverage, "\n\n\nCoverage (mean)",
                    color = color
                ),
                .barPlot(x, x$length, "\n\n\nLength", color = color),
                .barPlot(x, x$degeneracy, "\n\n\nDegeneracy", color = color)
            ),
            ncol = 6
        )
    }
}

.oligoFeaturePlot <- function(x) {
    if (all(x$type == "probe")) {
        patchwork::wrap_plots(list(
            .gcTmCoveragePlot(
                x[x$type == "probe", ],
                color = "grey20",
                type = "Probes"
            )
        ), ncol = 1)
    }
    else if (all(x$type == "primer")) {
        patchwork::wrap_plots(list(
            .gcTmCoveragePlot(
                x[x$type == "primer", ],
                color = "grey20",
                type = "Primers"
            )
        ), ncol = 1)
    } else {
        patchwork::wrap_plots(list(
            .gcTmCoveragePlot(
                x[x$type == "primer", ],
                color = "grey20", type = "Primers"
            ),
            .gcTmCoveragePlot(
                x[x$type == "probe", ],
                color = "grey20", type = "Probes"
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
            data = x[!duplicated(data.frame(x$start, x$end)), ],
            ggplot2::aes(
                xmin = start, xmax = end, ymin = 0.05, ymax = 0.65
            ), fill = "#424B54", alpha = 0.5
        ) +
        ggplot2::annotate(
            "label",
            x = start, y = 0.8,
            label = paste(
                "Assays n =", nrow(x)
            ), size = 3, hjust = 0, fontface = 2,
            color = "#424B54", fill = "white", label.size = NA
        ) +
        .themeRprimer(showYAxis = FALSE)
}

.assayFeaturePlot <- function(x, color = "grey20") {
    if (nrow(x) >= 10) {
        patchwork::wrap_plots(
            list(
                .boxPlot(
                    x, x$ampliconLength, "\n\nAmplicon length",
                    color = color
                ),
                .barPlot(
                    x, x$totalDegeneracy, "\n\nTotal degeneracy",
                    color = color
                )
            ),
            ncol = 2
        )
    } else {
        patchwork::wrap_plots(
            list(
                .barPlot(
                    x, x$ampliconLength, "\n\nAmplicon length",
                    color = color
                ),
                .barPlot(
                    x, x$totalDegeneracy, "\n\nTotal degeneracy",
                    color = color
                )
            ),
            ncol = 2
        )
    }
}

# Helpers for plotting an RprimerProfile =======================================

.identifyMask <- function(x) {
    x <- as.data.frame(x)
    x <- x[c("a", "c", "g", "t", "other", "gaps")]
    x <- as.matrix(x)
    masked <- apply(x, 1, function(y) all(is.na(y)))
    which(masked)
}

#' Calculate running average
#'
#' @param x A numeric vector.
#'
#' @param size
#' The number of observations in each average.
#' If \code{NULL}, the size will be set to the nearest positive, nonzero
#' integer to \code{length(x) / 100}.
#'
#' @return A data frame with position and running average of \code{x}.
#'
#' @noRd
.runningAverage <- function(x, size = NULL) {
    if (is.null(size)) {
        size <- round(length(x) / 100)
        if (size == 0) size <- 1
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

.identityPlot <- function(x, highlight = NULL, mask) {
    position <- identity <- average <- NULL
    averages <- .runningAverage(x$identity)
    xadj <- unique(x$position - seq_along(x$position))
    averages$position <- averages$position + xadj
    ggplot2::ggplot(
        data = x, ggplot2::aes(x = position, y = identity)
    ) +
        .highlightRegion(highlight) +
        ggplot2::geom_point(alpha = 1 / 3, shape = 1, color = "#93A8AC") +
        ggplot2::geom_line(
            data = averages, color = "#1B1C22",
            ggplot2::aes(x = position, y = average)
        ) +
        .maskRegion(mask + xadj) +
        ggplot2::ylim(0, 1) +
        ggplot2::ylab("Identity") +
        ggplot2::xlab("") +
        .themeRprimer(showXAxis = FALSE)
}

.entropyPlot <- function(x, highlight = NULL, mask) {
    position <- entropy <- average <- NULL
    averages <- .runningAverage(x$entropy)
    xadj <- unique(x$position - seq_along(x$position))
    averages$position <- averages$position + xadj
    ggplot2::ggplot(
        data = x, ggplot2::aes(x = position, y = entropy)
    ) +
        .highlightRegion(highlight) +
        ggplot2::geom_point(alpha = 1 / 3, shape = 1, color = "#93A8AC") +
        ggplot2::geom_line(
            data = averages,  color = "#1B1C22",
            ggplot2::aes(x = position, y = average)
        ) +
        .maskRegion(mask + xadj) +
        ggplot2::ylab("Entropy") +
        ggplot2::xlab("") +
        .themeRprimer(showXAxis = FALSE)
}

.gcPlot <- function(x, highlight = NULL, mask) {
    position <- average <- NULL
    gc <- ifelse(x$majority == "C" | x$majority == "G", 1, 0)
    gc[x$majority == "-"] <- 0.5
    gc[is.na(gc)] <- 0.5
    averages <- .runningAverage(gc)
    xadj <- unique(x$position - seq_along(x$position))
    averages$position <- averages$position + xadj
    ggplot2::ggplot(data = x, ggplot2::aes(x = position)) +
        ggplot2::geom_segment(
            color = "#93A8AC",
            ggplot2::aes(
                x = min(position), xend = max(position), y = 0.5, yend = 0.5
            )
        ) +
        .highlightRegion(highlight) +
        ggplot2::geom_line(
            data = averages, color = "#1B1C22",
            ggplot2::aes(x = position, y = average)
        ) +
        .maskRegion(mask + xadj) +
        ggplot2::xlab("") +
        ggplot2::ylab("GC-content") +
        ggplot2::ylim(0, 1) +
        .themeRprimer(showXAxis = FALSE)
}

.gapPlot <- function(x, highlight = NULL, mask) {
    position <- gaps <- NULL
    xadj <- unique(x$position - seq_along(x$position))
    ggplot2::ggplot(
        data = x, ggplot2::aes(x = position, y = gaps)
    ) +
        .highlightRegion(highlight) +
        ggplot2::geom_point(alpha = 1 / 3, shape = 1, color = "#93A8AC") +
        .maskRegion(mask + xadj) +
        ggplot2::ylim(0, 1) +
        ggplot2::xlab("Position") +
        ggplot2::ylab("Gaps") +
        .themeRprimer()
}

.maskRegion <- function(x) {
    ggplot2::geom_vline(xintercept = x, color = "grey80")
}

.highlightRegion <- function(highlight = NULL) {
    ggplot2::annotate(
        "rect",
        xmin = min(highlight), xmax = max(highlight), ymin = -Inf, ymax = Inf,
        color = "white", alpha = 0.4, fill = "#9B6A6C"
    )
}

.plotOverview <- function(x, highlight = NULL) {
    x <- as.data.frame(x)
    mask <- .identifyMask(x)
    patchwork::wrap_plots(
        list(
            .identityPlot(x, highlight, mask),
            .entropyPlot(x, highlight, mask),
            .gcPlot(x, highlight, mask),
            .gapPlot(x, highlight, mask)
        ),
        ncol = 1
    )
}

.plotNucleotides <- function(x, rc = FALSE) {
    position <- x$position
    x <- as.data.frame(x)
    select <- c("a", "c", "g", "t", "other")
    x <- x[names(x) %in% select]
    x <- t(as.matrix(x))
    colnames(x) <- position
    rownames(x)[-5] <- toupper(rownames(x))[-5]
    if (rc) {
        rownames(x) <- unname(lookup$complement[rownames(x)])
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
    if (rc) {
        x$Position <- factor(x$Position, levels = rev(levels(x$Position)))
    }
    basePalette <- c(
        "other" = "#D1D6DB",
        "A" = "#424B54", "C" = "#93A8AC", "G" = "#E2B4BD", "T" = "#9B6A6C"
    )
    ggplot2::ggplot(
        data = x, ggplot2::aes(x = Position, y = Frequency, fill = Base)
    ) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_manual(values = basePalette) +
        ggplot2::ylab("Proportion") +
        .themeRprimer(showLegend = TRUE)
}

# Helpers for plotting an RprimerMatchOligo object =============================

.plotMatch <- function(x, showLegend = TRUE, type = "oligo") {
    mismatches <- value <- NULL
    id <- as.factor(seq_len(nrow(x)))
    x <- x[!(grepl("id", names(x)))]
    x <- cbind(id, x)
    x <- suppressMessages({
        reshape2::melt(x)
    })
    names(x)[3] <- "mismatches"
    levels(x$mismatches) <- c(
        "0 mismatches", "1 mismatch",
        paste0(seq(2, length(levels(x$mismatches)) - 3), " mismatches"),
        paste0(length(levels(x$mismatches)) - 2, " or more mismatches"),
        "Off target matches"
    )
    if (type == "oligo") {
        yLabels <- x$iupacSequence
    } else {
        yLabels <- paste0(
            "assay ", id, ", ", type,
            " (length: ", nchar(x$iupacSequence), ")"
        )
    }
    onTarget <- x[x$mismatches != "Off target matches", ]
    offTarget <- x[x$mismatches == "Off target matches", ]
    ggplot2::ggplot(data = x, ggplot2::aes(x = id)) +
        ggplot2::geom_bar(
            data = onTarget,
            ggplot2::aes(x = id, y = value, fill = mismatches),
            stat = "identity", position = "stack"
        ) +
        ggplot2::geom_bar(
            data = offTarget,
            ggplot2::aes(x = id, y = value, fill = mismatches),
            stat = "identity", position = "stack", width = 0.1
        ) +
        ggplot2::ylab("Proportion") +
        ggplot2::xlab("") +
        ggplot2::scale_x_discrete(
            limits = rev(levels(factor(x$id))),
            labels = rev(yLabels)
        ) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(
            values = c(
                "#BDC9CC", "#B4B1B4", "#AC999C",
                "#A38183", "#9B6A6C", "#424B54"
            )
        ) +
        .themeRprimer(showLegend = showLegend)
}

# Helpers for plotting an RprimerMatchAssay object =============================

.plotMatchAssay <- function(x) {
    fwd <- x[grepl("Fwd$", names(x))]
    names(fwd) <- gsub("Fwd", "", names(fwd))
    rev <- x[grepl("Rev$", names(x))]
    names(rev) <- gsub("Rev", "", names(rev))
    if (any(grepl("Pr$", names(x)))) {
        pr <- x[grepl("Pr$", names(x))]
        names(pr) <- gsub("Pr", "", names(pr))
        patchwork::wrap_plots(
            list(
                .plotMatch(fwd, showLegend = FALSE, type = "fwd"),
                .plotMatch(pr, type = "pr", showLegend = FALSE),
                .plotMatch(rev, type = "rev")
            ),
            ncol = 3
        )
    } else {
        patchwork::wrap_plots(
            list(
                .plotMatch(fwd, showLegend = FALSE, type = "fwd"),
                .plotMatch(rev, type = "rev")
            ),
            ncol = 2
        )
    }
}

# Plot theme ===================================================================

.themeRprimer <- function(showXAxis = TRUE,
                          showYAxis = TRUE,
                          showLegend = FALSE) {
    if (showXAxis && showYAxis) {
        .showXYAxes(showLegend)
    } else if (showXAxis && !showYAxis) {
        .showXAxisHideYAxis(showLegend)
    } else if (!showXAxis && showYAxis) {
        .showYAxisHideXAxis(showLegend)
    } else {
        .hideXYAxes(showLegend)
    }
}

.showXYAxes <- function(showLegend = TRUE) {
    ggplot2::theme_bw() +
        ggplot2::theme(
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(size = 9),
            legend.position = ifelse(showLegend, "right", "none"),
            axis.title.y = ggplot2::element_text(
                size = 9,
                margin = ggplot2::margin(r = 10)
            ),
            axis.title.x = ggplot2::element_text(
                size = 9,
                margin = ggplot2::margin(t = 10)
            ),
            plot.title = ggplot2::element_text(size = 9),
            plot.margin = ggplot2::unit(rep(0.1, 4), "cm")
        )
}

.showXAxisHideYAxis <- function(showLegend = TRUE) {
    ggplot2::theme_bw() +
        ggplot2::theme(
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(size = 9),
            legend.position = ifelse(showLegend, "right", "none"),
            axis.title.x = ggplot2::element_text(
                size = 9,
                margin = ggplot2::margin(t = 10)
            ),
            panel.grid.major.y = ggplot2::element_blank(),
            panel.grid.minor.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(size = 9),
            plot.margin = ggplot2::unit(rep(0.1, 4), "cm")
        )
}

.showYAxisHideXAxis <- function(showLegend = TRUE) {
    ggplot2::theme_bw() +
        ggplot2::theme(
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(size = 9),
            legend.position = ifelse(showLegend, "right", "none"),
            axis.title.y = ggplot2::element_text(
                size = 9,
                margin = ggplot2::margin(r = 10)
            ),
            axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(size = 9),
            plot.margin = ggplot2::unit(rep(0.1, 4), "cm")
        )
}

.hideXYAxes <- function(showLegend = TRUE) {
    ggplot2::theme_bw() +
        ggplot2::theme(
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(size = 9),
            legend.position = ifelse(showLegend, "right", "none"),
            axis.title.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(size = 9),
            plot.margin = ggplot2::unit(rep(0.1, 4), "cm")
        )
}
