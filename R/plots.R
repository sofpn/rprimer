#' Plot the nucleotide distribution of an RprimerProfile object
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
#' plotNucleotides(exampleRprimerProfile[, 1:10]) ## Plot the first 10 bases
#' plotNucleotides(exampleRprimerProfile[, 1:10], rc = TRUE)
#' ## As reverse complement
#'
#' @export
setGeneric(
    "plotNucleotides", function(x, rc = FALSE) standardGeneric(
        "plotNucleotides"
    )
)

#' @describeIn plotNucleotides
#'
#' @export
setMethod("plotNucleotides", "RprimerProfile", function(x, rc = FALSE) {
    if (!(is.logical(rc))) {
        stop("'rc' must be set to 'TRUE' or 'FALSE'.", call. = FALSE)
    }
    x <- rpGetData(x)
    x <- x[rownames(x) != "-", ]
    if (rc == TRUE) {
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
        .themeRprimer(showXAxis = TRUE, showLegend = TRUE)
})

#' Plot an overview of an RprimerProfile object.
#'
#' @param x An \code{RprimerProfile} object.
#'
#' @return A plot.
#'
#' @export
setGeneric(
    "plotAlignmentOverview", function(x) standardGeneric(
        "plotAlignmentOverview"
    )
)

#' @describeIn plotAlignmentOverview
#'
#' @export
setMethod("plotAlignmentOverview", "RprimerProfile", function(x) {
    patchwork::wrap_plots(
        list(.identityPlot(x), .entropyPlot(x), .gcPlot(x), .gapPlot(x)),
        ncol = 1
    )
})

# Helpers =====================================================================

.identityPlot <- function(x) {
    position <- identity <- average <- NULL
    averages <- .runningAverage(x$identity)
    ggplot2::ggplot(
        data = x, ggplot2::aes(x = seq_along(position), y = identity)
    ) +
        ggplot2::geom_point(size = 0.3, ggplot2::aes(colour = identity)) +
        ggplot2::geom_line(
            data = averages,  color = "#404038",
            ggplot2::aes(x = position, y = average)
        ) +
        ggplot2::ylim(0, 1) +
        ggplot2::ylab("Identity") +
        ggplot2::xlab("") +
        ggplot2::scale_color_gradient(low = "#a9bac785", high = "#a9bac785") +
        .themeRprimer()
}

.entropyPlot <- function(x) {
    position <- entropy <- average <- NULL
    averages <- .runningAverage(x$entropy)
    ggplot2::ggplot(
        data = x, ggplot2::aes(x = seq_along(position), y = entropy)
    ) +
        ggplot2::geom_point(size = 0.3, ggplot2::aes(colour = entropy)) +
        ggplot2::geom_line(
            data = averages,  color = "#404038",
            ggplot2::aes(x = position, y = average)
        ) +
        ggplot2::ylab("Entropy") +
        ggplot2::xlab("") +
        ggplot2::scale_color_gradient(low = "#a9bac785", high = "#a9bac785") +
        .themeRprimer()
}

.gcPlot <- function(x) {
    position <- average <- NULL
    averages <- .gcRunningAverage(x$majority)
    ggplot2::ggplot(data = x, ggplot2::aes(x = seq_along(position))) +
        ggplot2::geom_segment(
            color = "#a9bac7",
            ggplot2::aes(x = 1, xend = length(position), y = 0.5, yend = 0.5)
        ) +
        ggplot2::geom_line(
            data = averages, color = "#404038",
            ggplot2::aes(
                x = position, y = average
            )
        ) +
        ggplot2::theme_light() +
        ggplot2::xlab("") +
        ggplot2::ylab("GC-content") +
        ggplot2::ylim(0, 1) +
        .themeRprimer()
}

.gapPlot <- function(x) {
    position <- gaps <- NULL
    xadj <- unique(x$position - seq_along(x$position))
    ggplot2::ggplot(
        data = x, ggplot2::aes(x = seq_along(position), y = gaps)
    ) +
        ggplot2::geom_point(size = 0.3, ggplot2::aes(color = gaps)) +
        ggplot2::ylim(0, 1) +
        ggplot2::xlab("Position") +
        ggplot2::scale_x_continuous(labels = function(x) x + xadj) +
        ggplot2::ylab("Gaps") +
        ggplot2::scale_color_gradient(low = "#a9bac785", high = "#a9bac785") +
        .themeRprimer(showXAxis = TRUE)
}

.themeRprimer <- function(showXAxis = FALSE, showLegend = FALSE) {
    if (showXAxis) {
        ggplot2::theme_light() +
            ggplot2::theme(
                legend.title = ggplot2::element_blank(),
                legend.position = ifelse(showLegend, "right", "none"),
                axis.title.y = ggplot2::element_text(
                    size = 9, color = "grey30",
                    margin = ggplot2::margin(r = 10)
                ),
                axis.title.x = ggplot2::element_text(
                    size = 9, color = "grey30",
                    margin = ggplot2::margin(t = 10)
                ),
                plot.margin = ggplot2::unit(rep(0, 4), "cm")
            )
    } else {
        ggplot2::theme_light() +
            ggplot2::theme(
                legend.title = ggplot2::element_blank(),
                legend.position = ifelse(showLegend, "right", "none"),
                axis.title.y = ggplot2::element_text(
                    size = 9, color = "grey30",
                    margin = ggplot2::margin(r = 10)
                ),
                axis.title.x = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                plot.margin = ggplot2::unit(rep(0, 4), "cm")
            )
    }
}

