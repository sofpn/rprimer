#' Plot an rprimer object (generic)
#'
#' @param x An rprimer object (see methods for details).
#'
#' @param ...
#' Additional arguments that should be passed to the plot, when
#' plotting an RprimerProfile object.
#'
#' @return A plot.
#'
#' @examples
#' data("exampleRprimerProfile")
#' rpPlot(exampleRprimerProfile[, 1:10]) ## Plot the first 10 bases
#' rpPlot(exampleRprimerProfile[, 1:10], rc = TRUE) ## As reverse complement
#'
#' data("exampleRprimerProperties")
#' rpPlot(exampleRprimerProperties)
#'
#' @export
setGeneric("rpPlot", function(x, ...) standardGeneric("rpPlot"))

#' @describeIn rpPlot Plot an RprimerProfile object.
#'
#' @param rc \code{TRUE/FALSE}. If the plotted sequence should be displayed
#' as reverse complement or not. The default is \code{FALSE}.
#'
#' @export
setMethod("rpPlot", "RprimerProfile", function(x, rc = FALSE) {
    if (!(is.logical(rc))) {
        stop("'rc' must be set to TRUE or FALSE", call. = FALSE)
    }
    x <- SummarizedExperiment::assay(x)
    x <- x[rownames(x) != "-", ]
    if (rc == TRUE) {
        rownames(x) <- unname(complementLookup[rownames(x)])
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
        ggplot2::theme_light() +
        ggplot2::theme(
            legend.title = ggplot2::element_blank()
        )
})

#' @describeIn rpPlot Plot an RprimerProperties object.
#'
#' @export
setMethod("rpPlot", "RprimerProperties", function(x) {
    x <- SummarizedExperiment::assay(x)
    cowplot::plot_grid(
        .identityPlot(x),
        .entropyPlot(x),
        .gcPlot(x),
        .gapPlot(x),
        align = "v",
        nrow = 4
    )
})

# Helpers =====================================================================

# Addera custom xlab har, dvs ej seq along, ska också va integer
# battre farger, transparent, tjockare linjer

#' Plot nucleotide identity
#'
#' @keywords internal
#'
#' @noRd
.identityPlot <- function(x) {
    Position <- Identity <- NULL
    averages <- .runningAverage(x$Identity)
    ggplot2::ggplot(
        data = x, ggplot2::aes(x = seq_along(Position), y = Identity)
        ) +
        ggplot2::geom_point(size = 0.2, ggplot2::aes(colour = Identity)) +
        ggplot2::geom_line(
            data = averages,
            ggplot2::aes(x = seq_along(Position), y = Average)
        ) +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_light() +
        ggplot2::xlab("") +
        ggplot2::theme(
            legend.position = "none",
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            plot.margin = ggplot2::unit(c(0, 1, 0, 0.5), "cm")
        )
}

#' Plot Shannon entropy
#'
#' @keywords internal
#'
#' @noRd
.entropyPlot <- function(x) {
    Position <- Entropy <- NULL
    averages <- .runningAverage(x$Entropy)
    ggplot2::ggplot(
        data = x, ggplot2::aes(x = seq_along(Position), y = Entropy)
        ) +
        ggplot2::geom_point(size = 0.2, ggplot2::aes(colour = Entropy)) +
        ggplot2::geom_line(
            data = averages,
            ggplot2::aes(x = seq_along(Position), y = Average)
        ) +
        ggplot2::theme_light() +
        ggplot2::ylim(0, NA) +
        ggplot2::xlab("") +
        ggplot2::theme(
            legend.position = "none",
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            plot.margin = ggplot2::unit(c(0, 1, 0, 0.5), "cm")
        )
}

#' Plot GC content
#'
#' @keywords internal
#'
#' @noRd
.gcPlot <- function(x) {
    Position <- NULL
    averages <- .gcRunningAverage(x$Majority)
    ggplot2::ggplot(data = x, ggplot2::aes(x = seq_along(Position))) +
        ggplot2::geom_line(
            data = averages,
            ggplot2::aes(x = seq_along(Position), y = Average)
        ) +
        ggplot2::geom_hline(yintercept = 0.5) +
        ggplot2::theme_light() +
        ggplot2::xlab("") +
        ggplot2::ylab("GC-content") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme(
            legend.position = "none",
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            plot.margin = ggplot2::unit(c(0, 1, 0, 0.5), "cm")
        )
}

#' Plot gap frequency
#'
#' @keywords internal
#'
#' @noRd
.gapPlot <- function(x) {
    Position <- Gaps <- NULL
    ggplot2::ggplot(data = x, ggplot2::aes(x = seq_along(Position), y = Gaps)) +
        ggplot2::geom_bar(stat = "identity", ggplot2::aes(color = Gaps)) +
        #ggplot2::scale_color_gradient(low = "blue", high = "red") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_light() +
        ggplot2::xlab("") +
        ggplot2::theme(
            legend.position = "none",
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            plot.margin = ggplot2::unit(c(0, 1, 0, 0.5), "cm")
        )
}
