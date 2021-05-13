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

.toDNAStringSetOligo <- function(x) {
    oligo <- x$sequence
    oligo <- unlist(lapply(seq_along(oligo), function(i) {
        names(oligo[[i]]) <- paste0(
            "oligo_", i, "_variant_", seq_along(oligo[[i]])
        ); oligo[[i]]
    }))
    Biostrings::DNAStringSet(oligo)
}

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

.toDNAStringSetAssay <- function(x) {
    fwd <- x$sequenceFwd
    fwd <- unlist(lapply(seq_along(fwd), function(i) {
        names(fwd[[i]]) <- paste0(
            "assay_", i, "_fwd_variant_", seq_along(fwd[[i]])
        ); fwd[[i]]
    }))
    fwd <- Biostrings::DNAStringSet(fwd)
    rev <- x$sequenceRev
    rev <- unlist(lapply(seq_along(rev), function(i) {
        names(rev[[i]]) <- paste0(
            "assay_", i, "_rev_variant_", seq_along(rev[[i]])
        ); rev[[i]]
    }))
    rev <- Biostrings::DNAStringSet(rev)
    rev <- Biostrings::reverseComplement(rev)
    if ("sequencePr" %in% names(x)) {
        pr <- x$sequencePr
        pr <- unlist(lapply(seq_along(pr), function(i) {
            names(pr[[i]]) <- paste0(
                "assay_", i, "_pr_variant_", seq_along(pr[[i]])
            ); pr[[i]]
        }))
        pr <- Biostrings::DNAStringSet(pr)
    } else {
        pr <- NULL
    }
    c(fwd, rev, pr)
}
