# Rprimer classes =============================================================

#' S4 classes for representation of different Rprimer objects
#'
#' @name Rprimer-classes
#'
#' @description
#' The rprimer package contains five different S4 classes. Each class is used
#' as input or output for the different functions within the oligo and assay
#' design workflow:
#'
#'  \itemize{
#'  \item \code{RprimerProfile}: output from \code{consensusProfile()},
#'  input for \code{oligos()}.
#'  \item \code{RprimerOligo}: output from \code{oligos()}, input
#'  for \code{assays()} and \code{checkMatch()}.
#'  \item \code{RprimerAssay}: output from \code{assays()}, input
#'  for \code{checkMatch()}.
#'  \item \code{RprimerMatchOligo}: output from \code{checkMatch()}.
#'  \item \code{RprimerMatchAssay}: output from \code{checkMatch()}.
#' }
#'
#' These classes extends the \code{S4Vectors::DataFrame} class
#' (Pages et al., 2020), without any additional slots.
#'
#' @section Coercion:
#' Each class can be converted to a traditional data frame, by using either
#' \code{as()} or \code{as.data.frame()}.
#'
#' Moreover, \code{as()} can also be used for converting oligo sequences
#' within an
#' \code{RprimerOligo} or \code{RprimerAssay} object into a
#' \code{Biostrings::DNAStringSet}
#' object (Pages et al., 2020).
#' Note that all oligo sequences will be written in the
#' same direction as the input alignment that was used to generate
#' the oligos.
#'
#' @param ...
#' A data frame or list to be converted into an Rprimer-object.
#'
#' @return
#' An Rprimer-object if validation succeeds, an error
#' message otherwise.
#'
#' @import methods
#'
#' @importClassesFrom S4Vectors DataFrame
#'
#' @importClassesFrom Biostrings DNAStringSet
#'
#' @seealso consensusProfile, oligos, assays, checkMatch
#'
#' @references
#' Pages, H., Lawrence, M., and Aboyoun, R. (2020). S4Vectors:
#' Foundation of vector-like and list-like containers in
#' Bioconductor. R package version 0.28.0.
#'
#' Pages, H., Aboyoun, P., Gentleman R., and DebRoy S. (2020). Biostrings:
#' Efficient manipulation of biological strings. R package version
#' 2.57.2.
#'
#' @examples
#' ## Constructors
#'
#' data("exampleRprimerProfile")
#' x <- as.data.frame(exampleRprimerProfile)
#' RprimerProfile(x)
#'
#' data("exampleRprimerOligo")
#' x <- as.data.frame(exampleRprimerOligo)
#' RprimerOligo(x)
#'
#' data("exampleRprimerAssay")
#' x <- as.data.frame(exampleRprimerAssay)
#' RprimerAssay(x)
#'
#' data("exampleRprimerMatchOligo")
#' x <- as.data.frame(exampleRprimerMatchOligo)
#' RprimerMatchOligo(x)
#'
#' data("exampleRprimerMatchAssay")
#' x <- as.data.frame(exampleRprimerMatchAssay)
#' RprimerMatchAssay(x)
#'
#' ## Coercion methods for RprimerOligo and RprimerAssay objects
#'
#' ## Convert an RprimerOligo object to a DNAStringSet
#' data("exampleRprimerOligo")
#'
#' ## Pick rows to convert
#' x <- exampleRprimerOligo[1:2, ]
#' as(x, "DNAStringSet")
#'
#' ## Convert an RprimerAssay object to a DNAStringSet
#' data("exampleRprimerAssay")
#'
#' ## Pick rows to convert
#' x <- exampleRprimerAssay[1:2, ]
#' as(x, "DNAStringSet")
NULL

# RprimerProfile ===============================================================

#' @rdname Rprimer-classes
.RprimerProfile <- setClass("RprimerProfile", contains = "DFrame")

#' @export
#'
#' @rdname Rprimer-classes
RprimerProfile <- function(...) {
    S4Vectors::DataFrame(..., row.names = NULL, check.names = TRUE) |>
        .RprimerProfile()
}

S4Vectors::setValidity2("RprimerProfile", \(object) {
    msg <- NULL
    colnames <- c(
        "position", "a", "c", "g", "t", "other", "gaps", "majority", "identity",
        "iupac", "coverage"
    )
    if (!all(colnames %in% names(object))) {
        msg <- c(
            msg, "This type of subsetting is not allowed for an
            RprimerProfile object."
        )
    }
    if (is.null(msg)) {
        TRUE
    } else {
        msg
    }
})

# RprimerOligo ================================================================

#' @rdname Rprimer-classes
.RprimerOligo <- setClass("RprimerOligo", contains = "DFrame")

#' @export
#'
#' @rdname Rprimer-classes
RprimerOligo <- function(...) {
    S4Vectors::DataFrame(..., row.names = NULL, check.names = TRUE) |>
        .RprimerOligo()
}

S4Vectors::setValidity2("RprimerOligo", \(object) {
    msg <- NULL
    colnames <- c(
        "type", "fwd", "rev", "start", "end", "length", "iupacSequence",
        "iupacSequenceRc", "identity",
        "coverage", "degeneracy", "gcContentMean", "gcContentRange", "tmMean",
        "tmRange", "deltaGMean", "deltaGRange",
        "sequence", "sequenceRc", "gcContent", "tm", "deltaG", "method",
        "score", "roiStart", "roiEnd"
    )
    if (!all(colnames %in% names(object))) {
        msg <- c(
            msg, "This type of subsetting is not allowed for an
            RprimerOligo object."
        )
    }
    if (is.null(msg)) {
        TRUE
    } else {
        msg
    }
})

# RprimerAssay ================================================================

#' @rdname Rprimer-classes
.RprimerAssay <- setClass("RprimerAssay", contains = "DFrame")

#' @export
#'
#' @rdname Rprimer-classes
RprimerAssay <- function(...) {
    S4Vectors::DataFrame(..., row.names = NULL, check.names = TRUE) |>
        .RprimerAssay()
}

S4Vectors::setValidity2("RprimerAssay", \(object) {
    msg <- NULL
    oligoColnames <- c(
        "start", "end",
        "length", "iupacSequence", "identity",
        "coverage", "degeneracy", "gcContentMean",
        "gcContentRange", "tmMean", "tmRange",
        "sequence",
        "gcContent", "tm", "deltaG",
        "method"
    )
    colnames <- c(
        "start", "end", "length",
        "totalDegeneracy", "score",
        paste0(oligoColnames, "Fwd"),
        paste0(oligoColnames, "Rev"),
        "roiStart", "roiEnd"
    )
    if (!all(colnames %in% names(object))) {
        msg <- c(
            msg, "This type of subsetting is not allowed for an
            RprimerAssay object."
        )
    }
    if (is.null(msg)) {
        TRUE
    } else {
        msg
    }
})

# RprimerMatchOligo ============================================================

#' @rdname Rprimer-classes
.RprimerMatchOligo <- setClass("RprimerMatchOligo", contains = "DFrame")

#' @export
#'
#' @rdname Rprimer-classes
RprimerMatchOligo <- function(...) {
    S4Vectors::DataFrame(..., row.names = NULL, check.names = TRUE) |>
        .RprimerMatchOligo()
}

S4Vectors::setValidity2("RprimerMatchOligo", \(object) {
    msg <- NULL
    colnames <- c(
        "perfectMatch", "oneMismatch", "twoMismatches",
        "threeMismatches", "fourOrMoreMismatches", "offTargetMatch"
    )
    if (!all(colnames %in% names(object))) {
        msg <- c(
            msg, "This type of subsetting is not allowed for an
            RprimerMatchOligo object."
        )
    }
    if (is.null(msg)) {
        TRUE
    } else {
        msg
    }
})

# RprimerMatchAssay ============================================================

#' @rdname Rprimer-classes
.RprimerMatchAssay <- setClass("RprimerMatchAssay", contains = "DFrame")

#' @export
#'
#' @rdname Rprimer-classes
RprimerMatchAssay <- function(...) {
    S4Vectors::DataFrame(..., row.names = NULL, check.names = TRUE) |>
        .RprimerMatchAssay()
}

S4Vectors::setValidity2("RprimerMatchAssay", \(object) {
    msg <- NULL
    oligoColnames <- c(
        "perfectMatch", "oneMismatch", "twoMismatches",
        "threeMismatches", "fourOrMoreMismatches", "offTargetMatch"
    )
    colnames <- c(
        paste0(oligoColnames, "Fwd"),
        paste0(oligoColnames, "Rev")
    )
    if (!all(colnames %in% names(object))) {
        msg <- c(
            msg, "This type of subsetting is not allowed for an
            RprimerMatchAssay object."
        )
    }
    if (is.null(msg)) {
        TRUE
    } else {
        msg
    }
})

# Coerce =======================================================================

#' @name coerce
#'
#' @rdname Rprimer-classes
setAs("RprimerOligo", "DNAStringSet", \(from) .toDNAStringSetOligo(from))

#' @name coerce
#'
#' @rdname Rprimer-classes
setAs("RprimerAssay", "DNAStringSet", \(from) .toDNAStringSetAssay(from))

# Helpers ======================================================================

.addNames <- function(x, type, additionalInfo = "") {
    unlist(lapply(seq_along(x), \(i) {
        names(x[[i]]) <- paste0(
            type, "_", i, additionalInfo, "_variant_", seq_along(x[[i]])
        )
        x[[i]]
    }))
}

.toDNAStringSetOligo <- function(x) {
    oligo <- x$sequence
    oligo <- .addNames(oligo, "oligo") |>
        Biostrings::DNAStringSet()
    oligo
}

.toDNAStringSetAssay <- function(x) {
    fwd <- x$sequenceFwd
    fwd <- .addNames(fwd, "assay", "_fwd") |>
        Biostrings::DNAStringSet()
    rev <- x$sequenceRev
    rev <- .addNames(rev, "assay", "_rev") |>
        Biostrings::DNAStringSet() |>
        Biostrings::reverseComplement()
    if ("sequencePr" %in% names(x)) {
        pr <- x$sequencePr
        pr <- .addNames(pr, "assay", "_pr") |>
            Biostrings::DNAStringSet()
    } else {
        pr <- NULL
    }
    c(fwd, rev, pr)
}
