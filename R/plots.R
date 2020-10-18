#' Plot nucleotide distribution
#'
#' @param x An \code{RprimerProfile} object.
#'
#' @param from
#' Optional. Start position (an integer).
#'
#' @param to
#' Optional. End position (an integer).
#'
#' @param rc
#' \code{TRUE} or {FALSE}. If the plotted sequence should be displayed
#' as reverse complement or not. Defaults to \code{FALSE}.
#'
#' @return A barplot with nucleotide distribution.
#'
#' @examples
#' data("exampleRprimerProfile")
#' plotNucleotides(exampleRprimerProfile, from = 1, to = 10)
#' ## Plot the first 10 bases
#' plotNucleotides(exampleRprimerProfile, from =  1, to = 10, rc = TRUE)
#' ## As reverse complement
#'
#' @export
plotNucleotides <-  function(x, from = NULL, to = NULL, rc = FALSE) {
 #   if (!methods::is(x, "RprimerProfile")) {
 #       stop("'x' must be an RprimerProfile object.", call. = FALSE)
 #   }
    if (!(is.logical(rc))) {
        stop("'rc' must be set to 'TRUE' or 'FALSE'.", call. = FALSE)
    }
    if (is.null(from)) from <- 1
    if (is.null(to)) to <- ncol(x)
    if (!is.numeric(from) || length(from) != 1) {
        stop("'from' must be a number.", call. = FALSE)
    }
    if (!is.numeric(to) || length(to) != 1) {
        stop("'to' must be a number.", call. = FALSE)
    }
    if (from >= to) {
        stop("'from' cannot be greater or equal to 'to'.", call. = FALSE)
    }
    x <- x[x$position >= from & x$position <= to, ]
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
        .themeRprimer(showXAxis = TRUE, showLegend = TRUE)
}

#' Plot sequence conservation, GC content and gap frequency
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
#' @export
plotConsensusProfile <- function(x, shadeFrom = NULL, shadeTo = NULL) {
   # if (!methods::is(x, "RprimerProfile")) {
#        stop("'x' must be an RprimerProfile object.", call. = FALSE)
#    }
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

# Helpers =====================================================================

.identityPlot <- function(x, shadeFrom = NULL, shadeTo = NULL) {
    position <- identity <- average <- NULL
    averages <- .runningAverage(x$identity)
    xadj <- unique(x$position - seq_along(x$position))
    shadeFrom <- shadeFrom - xadj
    shadeTo <- shadeTo - xadj
    ggplot2::ggplot(
        data = x, ggplot2::aes(x = seq_along(position), y = identity)
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
        .themeRprimer()
}

.entropyPlot <- function(x, shadeFrom = NULL, shadeTo = NULL) {
    position <- entropy <- average <- NULL
    averages <- .runningAverage(x$entropy)
    xadj <- unique(x$position - seq_along(x$position))
    shadeFrom <- shadeFrom - xadj
    shadeTo <- shadeTo - xadj
    ggplot2::ggplot(
        data = x, ggplot2::aes(x = seq_along(position), y = entropy)
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
        .themeRprimer()
}

.gcPlot <- function(x, shadeFrom = NULL, shadeTo = NULL) {
    position <- average <- NULL
    averages <- .gcRunningAverage(x$majority)
    xadj <- unique(x$position - seq_along(x$position))
    shadeFrom <- shadeFrom - xadj
    shadeTo <- shadeTo - xadj
    ggplot2::ggplot(data = x, ggplot2::aes(x = seq_along(position))) +
        ggplot2::geom_segment(
            color = "#a9bac7",
            ggplot2::aes(x = 1, xend = length(position), y = 0.5, yend = 0.5)
        ) +
        .shadeArea(shadeFrom = shadeFrom, shadeTo = shadeTo) +
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

.gapPlot <- function(x, shadeFrom = NULL, shadeTo = NULL) {
    position <- gaps <- NULL
    xadj <- unique(x$position - seq_along(x$position))
    shadeFrom <- shadeFrom - xadj
    shadeTo <- shadeTo - xadj
    ggplot2::ggplot(
        data = x, ggplot2::aes(x = seq_along(position), y = gaps)
    ) +
        .shadeArea(shadeFrom = shadeFrom, shadeTo = shadeTo) +
        ggplot2::geom_point(size = 0.3, ggplot2::aes(color = gaps)) +
        ggplot2::ylim(0, 1) +
        ggplot2::xlab("Position") +
        ggplot2::scale_x_continuous(labels = function(x) x + xadj) +
        ggplot2::ylab("Gaps") +
        ggplot2::scale_color_gradient(low = "#a9bac785", high = "#a9bac785") +
        .themeRprimer(showXAxis = TRUE)
}

.shadeArea <- function(shadeFrom = NULL, shadeTo = NULL) {
    ggplot2::annotate(
        "rect",
        xmin = shadeFrom, xmax = shadeTo, ymin = -Inf, ymax = Inf,
        color = "grey80", alpha = 0.2
    )
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
