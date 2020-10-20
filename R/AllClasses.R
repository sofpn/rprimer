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
#' \code{S4Vectors::DataFrame} constructor, without any additional arguments.
#'
#' @param ... A data frame or list to be converted into an \code{RprimerProfile}
#' object.
#'
#' @param row.names
#' \code{NULL} or a single integer or character string specifying a column to be
#' used as row names, or a character or integer vector giving the row names
#' for the data frame.
#'
#' @param check.names
#' logical. If \code{TRUE} then the names of the variables in the data frame are
#' checked to ensure that they are syntactically valid variable names and
#' are not duplicated.
#'
#' @param stringsAsFactors
#' logical. If character vectors should be converted to factors.
#'
#' @return An \code{RprimerAssay} object if validation succeeds, an error
#' message otherwise.
#'
#' @export
#'
#' @examples
#' data("exampleRprimerProfile")
#' x <- as.data.frame(exampleRprimerProfile)
#' RprimerProfile(x)
#' @importFrom S4Vectors DataFrame
RprimerProfile <- function(...,
                           row.names = NULL,
                           check.names = TRUE,
                           stringsAsFactors) {
    df <- DataFrame(...,
        row.names = row.names,
        check.names = check.names,
        stringsAsFactors = stringsAsFactors
    )
    .RprimerProfile(df)
}

# I did very simple validity checks here...
S4Vectors::setValidity2("RprimerProfile", function(object) {
    msg <- NULL
    #  if (names(object)[1] != "counts") {
    #    msg <- c(msg, "'counts' must be first assay")
    #  }
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
#' \code{S4Vectors::DataFrame} constructor, without any additional arguments.
#'
#' @describeIn RprimerOligo-class
#'
#' @param ... A data frame or list to be converted into an \code{RprimerProfile}
#' object.
#'
#' @param row.names
#' \code{NULL} or a single integer or character string specifying a column to be
#' used as row names, or a character or integer vector giving the row names
#' for the data frame.
#'
#' @param check.names
#' logical. If \code{TRUE} then the names of the variables in the data frame are
#' checked to ensure that they are syntactically valid variable names and
#' are not duplicated.
#'
#' @param stringsAsFactors
#' logical. If character vectors should be converted to factors.
#'
#' @return An \code{RprimerOligo} object if validation succeeds, an error
#' message otherwise.
#'
#' @export
#'
#' @examples
#' data("exampleRprimerOligo")
#' x <- as.data.frame(exampleRprimerOligo)
#' RprimerOligo(x)
#' @importFrom S4Vectors DataFrame
RprimerOligo <- function(...,
                         row.names = NULL,
                         check.names = TRUE,
                         stringsAsFactors) {
    df <- DataFrame(...,
        row.names = row.names,
        check.names = check.names,
        stringsAsFactors = stringsAsFactors
    )
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
#' \code{S4Vectors::DataFrame} constructor, without any additional arguments.
#'
#' @describeIn RprimerAssay-class
#'
#' @param ... A data frame or list to be converted into an \code{RprimerProfile}
#' object.
#'
#' @param row.names
#' \code{NULL} or a single integer or character string specifying a column to be
#' used as row names, or a character or integer vector giving the row names
#' for the data frame.
#'
#' @param check.names
#' logical. If \code{TRUE} then the names of the variables in the data frame are
#' checked to ensure that they are syntactically valid variable names and
#' are not duplicated.
#'
#' @param stringsAsFactors
#' logical. If character vectors should be converted to factors.
#'
#' @return An \code{RprimerAssay} object if validation succeeds, an error
#' message otherwise.
#'
#' @export
#'
#' @examples
#' data("exampleRprimerAssay")
#' x <- as.data.frame(exampleRprimerAssay)
#' RprimerAssay(x)
#' @importFrom S4Vectors DataFrame
RprimerAssay <- function(...,
                         row.names = NULL,
                         check.names = TRUE,
                         stringsAsFactors) {
    df <- DataFrame(...,
        row.names = row.names,
        check.names = check.names,
        stringsAsFactors = stringsAsFactors
    )
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
