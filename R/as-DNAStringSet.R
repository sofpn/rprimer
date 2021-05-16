#' Coerce an RprimerOligo object to a DNAStringset object
#'
#' @details
#' All sequences will be written in 5' to 3' direction.
#'
#' @describeIn RprimerOligo-class
#'
#' @export
#'
#' @import methods
setAs("RprimerOligo", "DNAStringSet", function(from) .toDNAStringSetOligo(from))

#' Coerce an RprimerAssay object to a DNAStringset object
#'
#' @details
#' All sequences will be written in 5' to 3' direction.
#'
#' @describeIn RprimerAssay-class
#'
#' @export
#'
#' @import methods
setAs("RprimerAssay", "DNAStringSet", function(from) .toDNAStringSetAssay(from))

# Helpers ======================================================================

.addNames <- function(x, type, additionalInfo = "") {
    unlist(lapply(seq_along(x), function(i) {
        names(x[[i]]) <- paste0(
            type, "_", i, additionalInfo, "_variant_", seq_along(x[[i]])
        )
        x[[i]]
    }))
}

.toDNAStringSetOligo <- function(x) {
    oligo <- x$sequence
    oligo <- .addNames(oligo, "oligo")
    Biostrings::DNAStringSet(oligo)
}

.toDNAStringSetAssay <- function(x) {
    fwd <- x$sequenceFwd
    fwd <- .addNames(fwd, "assay", "_fwd")
    fwd <- Biostrings::DNAStringSet(fwd)
    rev <- x$sequenceRev
    rev <- .addNames(rev, "assay", "_rev")
    rev <- Biostrings::DNAStringSet(rev)
    rev <- Biostrings::reverseComplement(rev)
    if ("sequencePr" %in% names(x)) {
        pr <- x$sequencePr
        pr <- .addNames(pr, "assay", "_pr")
        pr <- Biostrings::DNAStringSet(pr)
    } else {
        pr <- NULL
    }
    c(fwd, rev, pr)
}
