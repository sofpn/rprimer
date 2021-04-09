# tryCatch
# snyggare scoretabell i man,
# tester scores, tester match
# catch missing or partially missing targets!!!
# fix RprimerMatch class
# flexibilitet i plotten ocksa

#' Get indexes of perfectly matching and mismatching sequences
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerOligo")
#' data("exampleRprimerAlignment")
#' x <- exampleRprimerOligo$sequence[[1]]
#' .getMatchIndex(x, exampleRprimerAlignment)
.getMatchIndex <- function(x, target, maxMismatch = 3) {
    target <- Biostrings::DNAStringSet(target)
    x <- Biostrings::DNAStringSet(x)
    selection <- seq_along(target)
    res <- lapply(seq(0, maxMismatch), function(i) {
        result <- Biostrings::vcountPDict(x, target[selection], max.mismatch = i)
        selection <- which(colSums(result) == 0)
        selection
    })
    res[-length(res)] <- lapply(seq(1, length(res) - 1), function(i) {
        setdiff(res[[i]], res[[i + 1]])
    })
    res[[length(res) + 1]] <- setdiff(seq_along(target), unlist(res))
    res <- res[c(length(res), 1:(length(res) - 1))]
    res
}

#' Get number of perfectly matching and mismatching sequences
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerOligo")
#' data("exampleRprimerAlignment")
#' x <- head(exampleRprimerOligo)
#' .getMatchNumbers(x, exampleRprimerAlignment)
.getMatchNumbers <- function(x, target, maxMismatch = 3) {
    seq <- x$sequence
    numbers <- lapply(seq, function(x) {
        match <- .getMatchIndex(x, target, maxMismatch)
        match <- vapply(match, length, double(1L), USE.NAMES = FALSE)
    })
    numbers <- data.frame(do.call("rbind", numbers))
    names(numbers) <- c(
        paste0("n_", seq(0, maxMismatch), "_mm"),
        paste0("n_>=", maxMismatch + 1, "_mm")
    )
    cbind("iupacSequence" = x$iupacSequence, numbers)
}

.plotMatch <- function(x) {
    id <- as.character(seq_along(x$iupacSequence))
    x <- cbind(id, x)
    x <- reshape2::melt(x)
    names(x)[3] <- "mismatches"
    mismatchPalette <- grDevices::colorRampPalette(c("#c0cccf", "#435457"))
    levels(x$mismatches) <- c(
        "Perfect match", "1 mismatch",
        paste0(seq(2, length(levels(x$mismatches)) - 2), " mismatches"),
        paste0(length(levels(x$mismatches)) - 1, " or more mismatches")
    )
    ggplot2::ggplot(data = x, ggplot2::aes(
        fill = mismatches, x = id, y = value)
    ) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::ylab("Number of sequences") +
        ggplot2::xlab("") +
        ggplot2::scale_x_discrete(breaks = x$id, labels = x$iupacSequence) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(
            values = mismatchPalette(length(levels(x$mismatches)))
        ) +
        .themeRprimer(showLegend = TRUE)
}
