## Fix as IUPAC!

#' Convert to DNA string set (generic)
#'
#' \code{convertToDNAStringSet()} can be used for converting all sequence
#' variants of all oligos or assays to a \code{Biostrings::DNAStringSet} object.
#'
#' @param x
#' An \code{RprimerOligo} or \code{RprimerAssay} object.
#'
#' @param asRc
#' For \code{RprimerOligo} objects: if oligos should be
#' written as reverse complements.
#' For  \code{RprimerAssay} objects: if reverse primers and minus
#' sense probes should be
#' written as reverse complements.
#' Defaults to \code{TRUE}.
#'
#' @param asIUPAC
#' If sequences should be presented in IUPAC (degenerate) format. Defaults to
#' \code{TRUE}. If set to code{FALSE}, all sequence variants
#' of each oligo will be presented.
#'
#' @return A \code{Biostrings::DNAStringSet} object
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
#' convertToDNAStringSet(exampleRprimerAssay[1:2, ], asRc = FALSE)
setGeneric("convertToDNAStringSet", function(x, asRc = TRUE, asIUPAC = TRUE) standardGeneric("convertToDNAStringSet"))

#' Convert an RprimerOligo object to a DNAStringSet (method)
#'
#' @describeIn convertToDNAStringSet
#'
#' @export
setMethod("convertToDNAStringSet", "RprimerOligo", function(x,
                                                            asRc,
                                                            asIUPAC) {
    oligo <- if (asRc) x$sequenceRc else x$sequence
    oligo <- unlist(lapply(seq_along(oligo), function(i) {
        names(oligo[[i]]) <- paste0(
            "oligo_", i, "_variant_", seq_along(oligo[[i]])
        ); oligo[[i]]
    }))
    if (asRc) names(oligo) <- paste0(names(oligo), "_rc")
    Biostrings::DNAStringSet(oligo)
})

#' Convert an RprimerAssay object to a DNAStringSet (method)
#'
#' @describeIn convertToDNAStringSet
#'
#' @export
setMethod("convertToDNAStringSet", "RprimerAssay", function(x,
                                                            asRc,
                                                            asIUPAC) {
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
    if (asRc) {
        names(rev) <- paste0(names(rev), "_rc")
    }
    rev <- Biostrings::DNAStringSet(rev)
    if (!asRc) rev <- Biostrings::reverseComplement(rev)
    if ("sequencePr" %in% names(x)) {
        prPlus <- x$sequencePr
        prPlus <- unlist(lapply(seq_along(prPlus), function(i) {
            names(prPlus[[i]]) <- paste0(
                "assay_", i, "_pr_plus_variant_", seq_along(prPlus[[i]])
            ); prPlus[[i]]
        }))
        prPlus[!x$plusPr] <- "-"
        prPlus <- Biostrings::DNAStringSet(prPlus)
        if (asRc) {
            prMinus <- x$sequenceRcPr[x$minusPr]
        } else {
            prMinus <- x$sequencePr[x$minusPr]
        }
        prMinus <- unlist(lapply(seq_along(prMinus), function(i) {
            names(prMinus[[i]]) <- paste0(
                "assay_", i, "_pr_minus_variant_", seq_along(prMinus[[i]])
            ); prMinus[[i]]
        }))
        if (asRc) {
            names(prMinus) <- paste0(names(prMinus), "_rc")
        }
        prMinus[!x$minusPr] <- "-"
        prMinus <- Biostrings::DNAStringSet(prMinus)
    } else {
        prPlus <- prMinus <- NULL
    }
    out <- c(fwd, rev, prPlus, prMinus)
    out[out != "-"]
})
