## Kolla med iranges.............. "identify oligo target"
## tester
## assay

#' Check how oligos match to their target sequences (generic)
#'
#' @param x
#' An \code{RprimerOligo} or \code{RprimerAssay} object.
#'
#' @param target
#' The \code{Biostrings::DNAMultipleAlignment} object used for oligo/assay
#' design.
#'
#' ## Write details and limitations #########################
#' - only taerg region selected, seq w amb bases rmoved, no pos
#'
#'  the check will performed at all regions within the alignment except
#'  for the intended binding region (to detect any potential issues with off
#'  target binding).
#'  performed at the intended binding region (all other regions will be
#'  excluded to not generate false ). Here, to avoid false negative results,
#'  sequences with ambiguous bases at the oligo binding region are removed.
#'
#' The function is a wrapper to \code{Biostrings::vcountPDict()}
#' (Pages et al., 2020)
#'
#' @return
#' An \code{RprimerMatchOligo} or \code{RprimerMatchAssay} object. See details
#' below.
#'
#' @export
#'
#' @references
#' Pages, H., Aboyoun, P., Gentleman R., and DebRoy S. (2020). Biostrings:
#' Efficient manipulation of biological strings. R package version
#' 2.57.2.
setGeneric("checkMatch", function(x, target) standardGeneric("checkMatch"))

# Methods ======================================================================

#' Check match of an RprimerOligo object (method)
#'
#' @describeIn checkMatch
#'
#' @export
#'
#' @section Output (for oligos):
#'
#' The output contains the following information:
#'
#' \describe{
#'   \item{iupacSequence}{The oligo sequence in IUPAC format.}
#'   \item{perfectMatch}{Proportion of target sequences that matches perfectly
#'   to the oligo within the intended oligo binding region.}
#'   \item{oneMismatch}{Proportion of target sequences with one mismatch to
#'   the oligo
#'   within the intended oligo binding region.}
#'   \item{twoMismatches}{Proportion of target sequences with two mismatches
#'   to the oligo
#'   within the intended oligo binding region.}
#'   \item{threeMismatches}{Proportion of target sequences with
#'   three mismatches to the oligo
#'   within the intended oligo binding region.}
#'   \item{fourOrMoreMismatches}{Proportion of target sequences with
#'   four or more mismatches to the oligo
#'   within the intended oligo binding region.}
#'   \item{offTargetMatch}{Proportion of target sequences that matches
#'   to the oligo with
#'   no more than four mismatches to all other regions of the target alignment.}
#'  }
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

#' Check match of an RprimerAssay object (method)
#'
#' @describeIn checkMatch
#'
#' @export
#'
#' @examples
#' data("exampleRprimerAssay")
#' data("exampleRprimerAlignment")
#' x <- exampleRprimerAssay[1:10, ]
#' target <- exampleRprimerAlignment
#' checkMatch(x, target)
setMethod("checkMatch", "RprimerAssay", function(x, target) {
    if (!methods::is(target, "DNAMultipleAlignment")) {
        stop("'target' must be a DNAMultipleAlignment object.", call. = FALSE)
    }
    #match <- .checkMatchAssay(x, target)
    #RprimerMatchAssay(match)
})

# Helpers ======================================================================

.extractRange <- function(from, to, target, excludeAmbiguous = TRUE, ...) {
    selection <- target
    Biostrings::colmask(selection, ...) <- IRanges::IRanges(
        start = from, end = to
    )
    selection <- Biostrings::DNAStringSet(selection)
    if (excludeAmbiguous) {
        selection <- ShortRead::clean(selection)
    }
    selection
}


.getMatchIndex <- function(x, target) {
    res <- lapply(seq(0, 3), function(i) {
        result <- Biostrings::vcountPDict(
            x, target, max.mismatch = i
        )
        which(colSums(result) == 0)
    })
    res[-length(res)] <- lapply(seq(1, length(res) - 1), function(i) {
        setdiff(res[[i]], res[[i + 1]])
    })
    res[[length(res) + 1]] <- setdiff(seq_along(target), unlist(res))
    res[c(length(res), 1:(length(res) - 1))]
}

.getMatchIndexOffTarget <- function(x, target) {
    result <- Biostrings::vcountPDict(x, target, max.mismatch = 4)
    which(colSums(result) > 0)
}

.getMatchProportion <- function(x, target) {
    matching <- .getMatchIndex(x, target)
    matching <- vapply(matching, length, double(1L), USE.NAMES = FALSE)
    matching /  length(target)
}

.getMatchProportionOffTarget <- function(x, target) {
    matching <- .getMatchIndexOffTarget(x, target)
    length(matching) /  length(target)
}

.checkMatchOligos <- function(x, target) {
    onTarget <- lapply(seq_len(nrow(x)), function(i) {
        target <- .extractRange(
            x$start[[i]], x$end[[i]], target, invert = TRUE
        )
        check <- Biostrings::DNAStringSet(x$sequence[[i]])
        .getMatchProportion(check, target)
    })
    onTarget <- do.call("rbind", onTarget)
    onTarget <- as.data.frame(onTarget)

    offTarget <- lapply(seq_len(nrow(x)), function(i) {
        target <- .extractRange(
            x$start[[i]], x$end[[i]], target,
            excludeAmbiguous = FALSE, invert = FALSE
        )
        check <- Biostrings::DNAStringSet(x$sequence[[i]])
        .getMatchProportionOffTarget(check, target)
    })
    offTarget <- do.call("rbind", offTarget)
    offTarget <- as.data.frame(offTarget)
    all <- cbind( x$iupacSequence, onTarget, offTarget)
    names(all) <- c(
        "iupacSequence", "perfectMatch", "oneMismatch", "twoMismatches",
        "threeMismatches", "fourOrMoreMismatches", "offTargetMatch"
    )
    all
}
