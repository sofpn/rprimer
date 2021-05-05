# to do: write tests, example dataset and examples, validate class,
#  also make method for assay, add references
#######################################

#' Check how oligos match to a set of target sequences (generic)
#'
#' @param x
#' An \code{RprimerOligo} or \code{RprimerAssay} object.
#'
#' @param target
#' The \code{Biostrings::DNAMultipleAlignment} object used for generating the
#' consensus profile, with intended target sequences.
#'
#' ## Write details and limitations
#'
#' @return
#' A table with oligo sequences and proportion of sequences ######################
#' (An \code{RprimerMatchOligo} or \code{RprimerMatchAssay} object).
#'
#' @export
setGeneric("checkMatch", function(x, target) standardGeneric("checkMatch"))

# Methods ======================================================================

#' Check match of an RprimerOligo object (method)
#'
#' @describeIn checkMatch
#'
#' @export
#'
#' @examples
#' data("exampleRprimerOligo")
#' data("exampleRprimerAlignment")
#' x <- exampleRprimerOligo[1:10, ]
#' target <- exampleRprimerAlignment
#' checkMatch(x, target)
setMethod("checkMatch", "RprimerOligo", function(x, target) {
    if (!methods::is(target, "DNAMultipleAlignment")) {
        stop("'target' must be a DNAMultipleAlignment object.", call. = FALSE)
    }
    match <- .checkMatchOligos(x, target)
    RprimerMatchOligo(match)
})

# Helpers ======================================================================

.extractRange <- function(from, to, target) {
    selection <- target
    Biostrings::colmask(selection, invert = TRUE) <- IRanges::IRanges(
        start = from, end = to
    )
    selection <- Biostrings::DNAStringSet(selection)
    ShortRead::clean(selection) ## Remove sequences w wobble bases
}

.getMatchIndex <- function(x,
                           target,
                           maxMismatch = 3) {
    selection <- seq_along(target)
    res <- lapply(seq(0, maxMismatch), function(i) {
        result <- Biostrings::vcountPDict(
            x, target[selection], max.mismatch = i
        )
        which(colSums(result) == 0)
    })
    res[-length(res)] <- lapply(seq(1, length(res) - 1), function(i) {
        setdiff(res[[i]], res[[i + 1]])
    })
    res[[length(res) + 1]] <- setdiff(seq_along(target), unlist(res))
    res <- res[c(length(res), 1:(length(res) - 1))]
    names(res) <- c(
        paste0("n_", seq(0, maxMismatch), "_mm"),
        paste0("n_>=", maxMismatch + 1, "_mm")
    )
    res
}

.getMatchProportion <- function(x,
                                target,
                                maxMismatch = 3) {
    match <- .getMatchIndex(x, target, maxMismatch)
    numbers <- vapply(match, length, double(1L), USE.NAMES = FALSE)
    names(numbers) <- c(
        paste0("n", seq(0, maxMismatch), "mm"),
        paste0("n", maxMismatch + 1, "orMoreMm")
    )
   proportion <- numbers / length(target)
}

.checkMatchOligos <-function(x, target) {
    proportion <- lapply(seq_len(nrow(x)), function(i) {
        target <- .extractRange(x$start[[i]], x$end[[i]], target)
        check <- Biostrings::DNAStringSet(x$sequence[[i]])
        .getMatchProportion(check, target)
    })
    proportion <- do.call("rbind", proportion)
    proportion <- as.data.frame(proportion)
    cbind("iupacSequence" = x$iupacSequence, proportion)
}
