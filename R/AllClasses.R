# Profile - needs to be at least 10 bases or somethng
# positions must be consecutive

# oligo - type can only be primer or probe

# RprimerProfile ==============================================================

#' An S4 class for representation of a consensus profile
#'
#' @name RprimerProfile-class
#'
#' @description
#' \code{RprimerProfile} extends the \code{S4Vectors::DataFrame} class,
#' without any additional slots. It has the same accessors and
#' coercion, subsetting, and combining methods as the parent class, but has
#' some additional checks for validity.
#'
#' @export
#'
#' @import methods
#'
#' @importClassesFrom S4Vectors DataFrame
.RprimerProfile <- setClass("RprimerProfile", contains = "DataFrame")

#' RprimerProfile
#'
#' Constructs an \code{RprimerProfile} object in a similar fashion as the
#' \code{S4Vectors::DataFrame} constructor.
#'
#' @param ...
#' A data frame or list to be converted into an \code{RprimerProfile} object.
#'
#' @return
#' An \code{RprimerAssay} object if validation succeeds, an error
#' message otherwise.
#'
#' @export
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- as.data.frame(exampleRprimerProfile)
#' RprimerProfile(x)
#' @importFrom S4Vectors DataFrame
RprimerProfile <- function(...) {
    df <- DataFrame(..., row.names = NULL, check.names = TRUE)
    .RprimerProfile(df)
}

# to do..................
S4Vectors::setValidity2("RprimerProfile", function(object) {
    msg <- NULL
    colnames <- c(
        "position", "a", "c", "g", "t", "other", "gaps", "majority", "identity",
        "iupac", "entropy"
    )
    if (!all(colnames %in% names(object))) {
        msg <- c(
            msg, "The object must contain the following columns: \n
                 position, a, c, g, t, other, gaps, majority, identity."
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
#' without any additional slots. It has the same accessors and
#' coercion, subsetting, and combining methods as the parent class, but has
#' some additional checks for validity.
#'
#' @export
#'
#' @import methods
#'
#' @importClassesFrom S4Vectors DataFrame
.RprimerOligo <- setClass("RprimerOligo", contains = "DataFrame")

#' RprimerOligo
#'
#' Constructs an \code{RprimerOligo} object in a similar fashion as the
#' \code{S4Vectors::DataFrame} constructor.
#'
#' @describeIn RprimerOligo-class
#'
#' @param ...
#' A data frame or list to be converted into an \code{RprimerProfile} object.
#'
#' @return
#' An \code{RprimerOligo} object if validation succeeds, an error
#' message otherwise.
#'
#' @export
#'
#' @examples
#' data("exampleRprimerOligo")
#' x <- as.data.frame(exampleRprimerOligo)
#' RprimerOligo(x)
#' @importFrom S4Vectors DataFrame
RprimerOligo <- function(...) {
    df <- DataFrame(..., row.names = NULL, check.names = TRUE)
    .RprimerOligo(df)
}

S4Vectors::setValidity2("RprimerOligo", function(object) {
    msg <- NULL
    #  if (assayNames(object)[1] != "counts") {
    #    msg <- c(msg, "'counts' must be first assay")
    #  }
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
#' Constructs an \code{RprimerAssay} object in a similar fashion as the
#' \code{S4Vectors::DataFrame} constructor.
#'
#' @describeIn RprimerAssay-class
#'
#' @param ...
#' A data frame or list to be converted into an \code{RprimerProfile}
#' object.
#'
#' @return
#' An \code{RprimerAssay} object if validation succeeds, an error
#' message otherwise.
#'
#' @export
#'
#' @examples
#' data("exampleRprimerAssay")
#' x <- as.data.frame(exampleRprimerAssay)
#' RprimerAssay(x)
#' @importFrom S4Vectors DataFrame
RprimerAssay <- function(...) {
    df <- DataFrame(..., row.names = NULL, check.names = TRUE)
    .RprimerAssay(df)
}

S4Vectors::setValidity2("RprimerAssay", function(object) {
    msg <- NULL
    #  if (assayNames(object)[1] != "counts") {
    #    msg <- c(msg, "'counts' must be first assay")
    #  }
    if (is.null(msg)) {
        TRUE
    } else {
        msg
    }
})
