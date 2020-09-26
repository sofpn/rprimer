# Best practice:

# AllClasses.R
# AllGenerics.R

# allow subsetting but names must be kept! Also drop = F!
# getAlignmentProfile

#  New class ==================================================================

#' @export
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.RprimerProfile <- setClass(
  "RprimerProfile", contains = "SummarizedExperiment"
)

# Constructor =================================================================

#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
RprimerProfile <- function(x, ...) {
  se <- SummarizedExperiment(list(x = x), ...)
  .RprimerProfile(se)
}

# Validity method =============================================================

S4Vectors::setValidity2("RprimerProfile", function(object) {
  msg <- NULL
  if (!is.matrix(SummarizedExperiment::assay(object))) {
    msg <- c(msg, "The object is not a matrix.")
  }
  if (!is.numeric(SummarizedExperiment::assay(object))) {
    msg <- c(msg, "The object is not in numeric format.")
  }
 if (min(SummarizedExperiment::assay(object)) < 0) {
    msg <- c(msg, "The object contains values < 0.")
  }
  if (max(SummarizedExperiment::assay(object)) > 1) {
    msg <- c(msg, "The object contains values > 1.")
  }
  if (is.null(rownames(SummarizedExperiment::assay(object)))) {
    msg <- c(msg, "The object does not have rownames.")
  }
  #if (any(grepl(paste0("[^", dnaBases, "]"), row.names(SummarizedExperiment::assay(object))))) {
  #  msg <- c(msg, "The object has invalid rownames.")
  #}
  if (is.null(msg)) {
    TRUE
  } else {
    msg
  }
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
