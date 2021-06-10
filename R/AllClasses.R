# RprimerProfile ===============================================================

#' An S4 class for representation of a consensus profile
#'
#' @name RprimerProfile-class
#'
#' @description
#' \code{RprimerProfile} extends the \code{S4Vectors::DataFrame} class,
#' without any additional slots.
#'
#' @export
#'
#' @import methods
#'
#' @importClassesFrom S4Vectors DataFrame
.RprimerProfile <- setClass("RprimerProfile", contains = "DataFrame")

#' RprimerProfile
#'
#' The constructor, \code{RprimerProfile()},
#' constructs an \code{RprimerProfile} object in a similar fashion as the
#' \code{S4Vectors::DataFrame()} constructor.
#'
#' @param ...
#' A data frame or list to be converted into an \code{RprimerProfile} object.
#'
#' @return
#' An \code{RprimerProfile} object if validation succeeds, an error
#' message otherwise.
#'
#' @describeIn RprimerProfile-class
#'
#' @export
#'
#' @importFrom S4Vectors DataFrame
#'
#' @references
#' Pages, H., Lawrence, M., and Aboyoun, R. (2020). S4Vectors:
#' Foundation of vector-like and list-like containers in
#' Bioconductor. R package version 0.28.0.
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- as.data.frame(exampleRprimerProfile)
#' RprimerProfile(x)
RprimerProfile <- function(...) {
    df <- DataFrame(..., row.names = NULL, check.names = TRUE)
    .RprimerProfile(df)
}

S4Vectors::setValidity2("RprimerProfile", function(object) {
    msg <- NULL
    colnames <- c(
        "position", "a", "c", "g", "t", "other", "gaps", "majority", "identity",
        "iupac", "entropy", "coverage"
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

#' An S4 class for representation of oligos
#'
#' @name RprimerOligo-class
#'
#' @description
#' \code{RprimerOligo} extends the \code{S4Vectors::DataFrame} class,
#' without any additional slots.
#'
#' @export
#'
#' @import methods
#'
#' @importClassesFrom S4Vectors DataFrame
.RprimerOligo <- setClass("RprimerOligo", contains = "DataFrame")

#' RprimerOligo
#'
#' The constructor, \code{RprimerOligo()},
#' constructs an \code{RprimerOligo} object in a similar fashion as the
#' \code{S4Vectors::DataFrame()} constructor.
#'
#' @describeIn RprimerOligo-class
#'
#' @param ...
#' A data frame or list to be converted into an \code{RprimerOligo} object.
#'
#' @return
#' An \code{RprimerOligo} object if validation succeeds, an error
#' message otherwise.
#'
#' @export
#'
#' @importFrom S4Vectors DataFrame
#'
#' @references
#' Pages, H., Lawrence, M., and Aboyoun, R. (2020). S4Vectors:
#' Foundation of vector-like and list-like containers in
#' Bioconductor. R package version 0.28.0.
#'
#' @examples
#' data("exampleRprimerOligo")
#' x <- as.data.frame(exampleRprimerOligo)
#' RprimerOligo(x)
RprimerOligo <- function(...) {
    df <- DataFrame(..., row.names = NULL, check.names = TRUE)
    .RprimerOligo(df)
}

S4Vectors::setValidity2("RprimerOligo", function(object) {
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

#' An S4 class for representation of PCR assays
#'
#' @name RprimerAssay-class
#'
#' @description
#' \code{RprimerAssay} extends the \code{S4Vectors::DataFrame} class,
#' without any additional slots. It has the same accessors and
#' coercion, subsetting, and combining methods as the parent class, but has
#' some additional checks for validity.
#'
#' @export
#'
#' @import methods
#'
#' @importClassesFrom S4Vectors DataFrame
.RprimerAssay <- setClass("RprimerAssay", contains = "DataFrame")

#' RprimerAssay
#'
#' The constructor, \code{RprimerAssay()},
#' constructs an \code{RprimerAssay} object in a similar fashion as the
#' \code{S4Vectors::DataFrame()} constructor.
#'
#' @describeIn RprimerAssay-class
#'
#' @param ...
#' A data frame or list to be converted into an \code{RprimerAssay}
#' object.
#'
#' @return
#' An \code{RprimerAssay} object if validation succeeds, an error
#' message otherwise.
#'
#' @export
#'
#' @importFrom S4Vectors DataFrame
#'
#' @references
#' Pages, H., Lawrence, M., and Aboyoun, R. (2020). S4Vectors:
#' Foundation of vector-like and list-like containers in
#' Bioconductor. R package version 0.28.0.
#'
#' @examples
#' data("exampleRprimerAssay")
#' x <- as.data.frame(exampleRprimerAssay)
#' RprimerAssay(x)
RprimerAssay <- function(...) {
    df <- DataFrame(..., row.names = NULL, check.names = TRUE)
    .RprimerAssay(df)
}

S4Vectors::setValidity2("RprimerAssay", function(object) {
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

#' An S4 class for representation of match proportions of oligos
#'
#' @name RprimerMatchOligo-class
#'
#' @description
#' \code{RprimerMatchOligo} extends the \code{S4Vectors::DataFrame} class,
#' without any additional slots. It has the same accessors and
#' coercion, subsetting, and combining methods as the parent class, but has
#' some additional checks for validity.
#'
#' @export
#'
#' @import methods
#'
#' @importClassesFrom S4Vectors DataFrame
.RprimerMatchOligo <- setClass("RprimerMatchOligo", contains = "DataFrame")

#' RprimerMatchOligo
#'
#' The constructor, \code{RprimerMatchOligo()},
#' constructs an \code{RprimerMatchOligo} object in a similar fashion as the
#' \code{S4Vectors::DataFrame()} constructor.
#'
#' @describeIn RprimerMatchOligo-class
#'
#' @param ...
#' A data frame or list to be converted into an \code{RprimerMatchOligo}
#' object.
#'
#' @return
#' An \code{RprimerMatchOligo} object if validation succeeds, an error
#' message otherwise.
#'
#' @export
#'
#' @importFrom S4Vectors DataFrame
#'
#' @references
#' Pages, H., Lawrence, M., and Aboyoun, R. (2020). S4Vectors:
#' Foundation of vector-like and list-like containers in
#' Bioconductor. R package version 0.28.0.
#'
#' @examples
#' data("exampleRprimerMatchOligo")
#' x <- as.data.frame(exampleRprimerMatchOligo)
#' RprimerMatchOligo(x)
RprimerMatchOligo <- function(...) {
    df <- DataFrame(..., row.names = NULL, check.names = TRUE)
    .RprimerMatchOligo(df)
}

S4Vectors::setValidity2("RprimerMatchOligo", function(object) {
    msg <- NULL
    colnames <- c(
        "perfectMatch", "oneMismatch", "twoMismatches",
        "threeMismatches", "fourOrMoreMismatches"
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

#' An S4 class for representation of match percentages for assays
#'
#' @name RprimerMatchAssay-class
#'
#' @description
#' \code{RprimerMatchAssay} extends the \code{S4Vectors::DataFrame} class,
#' without any additional slots. It has the same accessors and
#' coercion, subsetting, and combining methods as the parent class, but has
#' some additional checks for validity.
#'
#' @export
#'
#' @import methods
#'
#' @importClassesFrom S4Vectors DataFrame
.RprimerMatchAssay <- setClass("RprimerMatchAssay", contains = "DataFrame")

#' RprimerMatchAssay
#'
#' The constructor, \code{RprimerMatchAssay()},
#' constructs an \code{RprimerMatchAssay} object in a similar fashion as the
#' \code{S4Vectors::DataFrame()} constructor.
#'
#' @describeIn RprimerMatchAssay-class
#'
#' @param ...
#' A data frame or list to be converted into an \code{RprimerMatchAssay}
#' object.
#'
#' @return
#' An \code{RprimerMatchAssay} object if validation succeeds, an error
#' message otherwise.
#'
#' @export
#'
#' @importFrom S4Vectors DataFrame
#'
#' @references
#' Pages, H., Lawrence, M., and Aboyoun, R. (2020). S4Vectors:
#' Foundation of vector-like and list-like containers in
#' Bioconductor. R package version 0.28.0.
#'
#' @examples
#' data("exampleRprimerMatchAssay")
#' x <- as.data.frame(exampleRprimerMatchAssay)
#' RprimerMatchAssay(x)
RprimerMatchAssay <- function(...) {
    df <- DataFrame(..., row.names = NULL, check.names = TRUE)
    .RprimerMatchAssay(df)
}

S4Vectors::setValidity2("RprimerMatchAssay", function(object) {
    msg <- NULL
    oligoColnames <- c(
        "perfectMatch", "oneMismatch", "twoMismatches",
        "threeMismatches", "fourOrMoreMismatches"
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

#' Coerce an RprimerOligo or RprimerAssay object to a DNAStringSet object
#'
#' \code{as} can be used for converting oligo sequences within an
#' RprimerOligo or RprimerAssay object into a DNAStringSet object
#' (Pages et al., 2020).
#'
#' @name coerce
#'
#' @aliases coerce, RprimerOligo, RprimerAssay
#'
#' @export
#'
#' @import methods
#'
#' @importClassesFrom Biostrings DNAStringSet
#'
#' @references
#' Pages, H., Aboyoun, P., Gentleman R., and DebRoy S. (2020). Biostrings:
#' Efficient manipulation of biological strings. R package version
#' 2.57.2.
#'
#' @examples
#' ## Convert an RprimerOligo object to a DNAStringSet
#' data("exampleRPrimerOligo")
#'
#' ## Pick rows to convert
#' x <- exampleRprimerOligo[1:2, ]
#' as(x, DNAStringSet)
#'
#' ## Convert an RprimerAssay object to a DNAStringSet
#' data("exampleRPrimerAssay")
#'
#' ## Pick rows to convert
#' x <- exampleRprimerAssay[1:2, ]
#' as(x, DNAStringSet)
setAs("RprimerOligo", "DNAStringSet", function(from) .toDNAStringSetOligo(from))

#'  @describeIn coerce
#'
#'  @export
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
