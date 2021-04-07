# tryCatch
# snyggare scoretabell i man,
# tester scores, tester match

# catch missing or partially missing targets!!! ......................
# fix RprimerMatch class

# "kor slut" pa antalet mismatches
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
.getMatchIndex <- function(x, target, maxMismatch = 6) {
    target <- Biostrings::DNAStringSet(target)
    x <- Biostrings::DNAStringSet(x)
    res <- lapply(seq(0, maxMismatch), function(i) {
        result <- Biostrings::vcountPDict(x, target, max.mismatch = i)
        which(colSums(result) == 0)
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
#' .getMatchProportion(x, exampleRprimerAlignment)
.getMatchProportion <- function(x, target, maxMismatch = 3) {
    seq <- x$sequence
    proportions <- lapply(seq, function(x) {
        match <- .getMatchIndex(x, target, maxMismatch)
        match <- vapply(match, function(x) {
            length(x) / length(target)
        }, double(1L), USE.NAMES = FALSE)
    })
    proportions <- data.frame(do.call("rbind", proportions))
    cbind("iupacSequence" = x$iupacSequence, proportions)
}

.plotMatch <- function(x) {
    x <- reshape2::melt(x)
    names(x)[2] <- "mismatches"
    ggplot2::ggplot(data = x, ggplot2::aes(
        fill = mismatches, x = iupacSequence, y = value)
    ) +
        ggplot2::geom_bar(stat = "identity", position = "fill") +
        #ggplot2::scale_fill_manual("Perfect match" = "green") +
        ggplot2::coord_flip() +
        scale_fill_discrete(
            name = "mismatch", labels = c(
                "Perfect match", "1 mismatch", "2 mismatches", "3 mismatches",
                "4 or more mismatches")
            ) +
        ggplot2::theme(
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            legend.title = ggplot2::element_blank()
        )
}


