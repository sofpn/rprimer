#' Plot an Rprimer object (generic)
#'
#' @param x An Rprimer object (see methods for details).
#'
#' @param ...
#' Additional arguments that should be passed to the plot, when
#' plotting an \code{RprimerProfile} object.
#'
#' @return A plot.
#'
#' @examples
#'
#' data("exampleRprimerProfile")
#' rpPlot(exampleRprimerProfile[, 1:10]) ## Plot the first 10 bases
#' rpPlot(exampleRprimerProfile[, 1:10], rc = TRUE) ## As reverse complement
#'
#' @export
setGeneric("rpPlot", function(x, ...) standardGeneric("rpPlot"))

#' @describeIn rpPlot Plot an \code{RprimerProfile} object.
#'
#' @param rc \code{TRUE/FALSE}. If the plotted sequence should be displayed
#' as reverse complement or not. The default is \code{FALSE}.
#'
#' @export
setMethod("rpPlot", "RprimerProfile", function(x, rc = FALSE) {
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
        ggplot2::theme_light() +
        ggplot2::theme(
            legend.title = ggplot2::element_blank()
        )
})
