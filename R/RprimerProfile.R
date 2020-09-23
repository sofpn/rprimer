# Best practice:

# AllClasses.R
# AllGenerics.R

# allow subsetting by col but not row

#  New class ==================================================================

#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.RprimerProfile <- setClass(
  "RprimerProfile", contains = "SummarizedExperiment"
)

# Constructor =================================================================

#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
RprimerProfile <- function(x, ...) {
  x <- x[, colSums(!is.na(x)) > 0]
  x <- x[(rownames(x) != "+" & rownames(x) != "."), ]
  se <- SummarizedExperiment(list(x = x), ...)
  .RprimerProfile(se)
}

# Validity method =============================================================

S4Vectors::setValidity2("RprimerProfile", function(object) {
  msg <- NULL
  if (!is.matrix(assay(object))) {
    msg <- c(msg, "The object is not a matrix.")
  }
  if (!is.numeric(assay(object))) {
    msg <- c(msg, "The object is not in numeric format.")
  }
 if (min(assay(object)) < 0) {
    msg <- c(msg, "The object contains values < 0.")
  }
  if (max(assay(object)) > 1) {
    msg <- c(msg, "The object contains values > 1.")
  }
  if (is.null(rownames(assay(object)))) {
    msg <- c(msg, "The object does not have rownames.")
  }
  if (nrow(assay(object)) != 16) {
    msg <- c("The object does not have 16 rows.")
  }
  if (any(grepl(paste0("[^", allBases, "]"), row.names(assay(object))))) {
    msg <- c(msg, "The object has invalid rownames.")
  }
  if (is.null(msg)) {
    TRUE
  } else {
    msg
  }
})

# Getter method ===============================================================

#' @export
setGeneric("rpShow", function(x, ...) standardGeneric("rpShow"))

#' @export
#' @importFrom SummarizedExperiment
setMethod("rpShow", "RprimerProfile", function(x, withDimNames = TRUE) {
  assay(x, withDimnames = withDimNames)
})

# Generics and methods ========================================================

#' @export
setGeneric("rpPlot", function(x, ...) standardGeneric("rpPlot"))


#' @export
setMethod("rpPlot", "RprimerProfile", function(x, rc = FALSE) {
  rc = rc

  plot(x)
})

# Inherits ====================================================================

#' Check if an object is as RprimerProfile-object
#'
#' @param x An 'RprimerProfile'-like object.
#'
#' @return \code{TRUE} or \code{FALSE}.
#'
#' @export
is.RprimerProfile <- function(x) inherits(x, "RprimerProfile")
