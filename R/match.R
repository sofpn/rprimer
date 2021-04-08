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
#' x <- head(exampleRprimerOligo)
#' .getMatchIndex(x, exampleRprimerAlignment)
.getMatchIndex <- function(x, target, maxMismatch = 3) {
    target <- Biostrings::DNAStringSet(target)
    x <- Biostrings::DNAStringSet(x)
    res <- lapply(seq(0, maxMismatch), function(i) {
        result <- Biostrings::vcountPDict(x, target, max.mismatch = i)
        result <- which(colSums(result) == 0)
    })
    res[seq_len(maxMismatch - 1)] <- lapply(
        seq_len(maxMismatch - 1), function(i) {
        res[[i]] <- setdiff(res[[i]], res[[i + 1]])
        }
    )
    names(res) <- c(
        paste0("n_", seq_len(maxMismatch), "_mm"),
        paste0("n_", maxMismatch, "or_more_mm")
    )
    res$pm <- setdiff(seq_along(target), unlist(res))
    res <- res[c(length(res), 1:(length(res) - 1))]
    res
}

#' Get proportion of perfectly matching and mismatching sequences
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
    cbind("iupacSequence" = x$iupacSequence, numbers)
}

.plotMatch <- function(x) {
    id <- as.character(seq_len(nrow(x)))
    x <- cbind(x, id)
    x <- reshape2::melt(x)
    names(x)[2] <- "mismatches"
    levels(x$mismatches) <- c(
        "Perfect match", "1 mismatch", "2 mismatches", "3 mimsathces",
        "4 or more mismatches"
    )
    mismatchPalette <- rev(c(
        "#435457", "#586f73", "#6d898f", "#89a0a4", "#c0cccf"
    ))
    ggplot2::ggplot(data = x, ggplot2::aes(
        fill = mismatches, x = id, y = value)
    ) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::ylab("Number of sequences") +
        ggplot2::xlab("") +
       # ggplot2::scale_x_discrete(x$iupacSequence) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = mismatchPalette) +
        .themeRprimer(showLegend = TRUE)
}





