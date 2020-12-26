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
        stop("'rc' must be set to TRUE or FALSE.", call. = FALSE)
    }
    position <- x$position
    x <- as.data.frame(x)
    x <- dplyr::select(x, c("a", "c", "g", "t", "other"))
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
        "other" = "#E8E9ED", "A" = "#FFC759", "C" = "#FF7B9C",
        "G" = "#607196", "T" = "#BABFD1"
    )
    ggplot2::ggplot(
        data = x, ggplot2::aes(x = Position, y = Frequency, fill = Base)
    ) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_manual(values = basePalette) +
        .themeRprimer(showLegend = TRUE, rotateX = TRUE)
}
