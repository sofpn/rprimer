#' Convert to DNA string set (generic)
#'
#' \code{convertToDNAStringSet()} can be used for converting all sequence
#' variants of all oligos or assays to a \code{Biostrings::DNAStringSet} object.
#'
#' @param x
#' An \code{RprimerProfile}, \code{RprimerOligo} or \code{RprimerAssay} object.
#'
#' @param revAsRc
#' If reverse primers/negative sense probes should be written as reverse
#' complements. Defaults to \code{TRUE}.
#'
#' @return A Biostrings::DNAStringSet
#'
#' @export
#'
#' @references
#' H. Pages, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings:
#' Efficient manipulation of biological strings. R package version
#' 2.57.2.
#'
#' @examples
#' data("exampleRprimerOligo")
#' convertToDNAStringSet(exampleRprimerOligo[1:2, ])
#'
#' data("exampleRprimerAssay")
#' convertToDNAStringSet(exampleRprimerAssay[1:2, ], revAsRc = FALSE)
setGeneric("convertToDNAStringSet", function(x, revAsRc = TRUE) standardGeneric("convertToDNAStringSet"))

#' Convert an RprimerOligo object to a DNAStringSet (method)
#'
#' @describeIn convertToDNAStringSet
#'
#' @export
setMethod("convertToDNAStringSet", "RprimerOligo", function(x,
                                                            revAsRc) {
    fwd <- x$sequence[x$fwd]
    fwd <- unlist(lapply(seq_along(fwd), function(i) {
        names(fwd[[i]]) <- paste0(
            "fwd_", i, "_variant_", seq_along(fwd[[i]])
        ); fwd[[i]]
    }))
    rev <- if (revAsRc) x$sequenceRc[x$rev] else x$sequence[x$rev]
    rev <- unlist(lapply(seq_along(rev), function(i) {
        names(rev[[i]]) <- paste0(
            "rev_", i, "_variant_", seq_along(rev[[i]])
        ); rev[[i]]
    }))
    Biostrings::DNAStringSet(c(fwd, rev))
})

#' Convert an RprimerAssay object to a DNAStringSet (method)
#'
#' @describeIn convertToDNAStringSet
#'
#' @export
setMethod("convertToDNAStringSet", "RprimerAssay", function(x,
                                                            revAsRc) {
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
    if (!revAsRc) rev <- Biostrings::reverseComplement(rev)
    if ("sequencePr" %in% names(x)) {
        prPlus <- x$sequencePr
        prPlus <- unlist(lapply(seq_along(pr), function(i) {
            names(pr[[i]]) <- paste0(
                "assay_", i, "_pr_variant_", seq_along(pr[[i]]) ########### IF pr minus!!!!!!!!!!!!!
            ); pr[[i]]
        }))
        prMinus <- Biostrings::DNAStringSet(prPlus)
        prMinus <- x$sequenceRcPr
        prMinus <- unlist(lapply(seq_along(prMinus), function(i) {
            names(prMinus[[i]]) <- paste0(
                "assay_", i, "_pr_rc_variant_", seq_along(prMinus[[i]])
            ); prMinus[[i]]
        }))
        prMinus <- Biostrings::DNAStringSet(prMinus)
    } else {
        prPlus <- prMinus <- NULL
    }
    c(fwd, rev, prPlus, prMinus)
})
