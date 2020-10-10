#' An S4 class for representation of an alignment profile
#'
#' An extension of the S4 class SummarizedExperiment::SummarizedExperiment.
#'
#' A numeric matrix with the proportion of each
#' nucleotide at each position within an alignment
#' of DNA sequences. The matrix has six rows,
#' named 'A', 'C', 'G', 'T', '-' and 'Other'. '-' represents gaps and
#' 'Other' represents nucleotides other than A, C, G and T, e.g.
#' wobble bases. The columns are named
#' according to which position they correspond to in the alignment.
#'
#' @name RprimerProfile-class
#'
#' @rdname RprimerProfile-class
#'
#' @export
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setClass("RprimerProfile", contains = "SummarizedExperiment")

# Constructor =================================================================

#' @describeIn RprimerProfile-class
#'
#' @param x An object.
#'
#' @param ... Additional arguments.
#'
#' @return An RprimerProfile object if validation succeeds. An error
#' message if not.
#'
#' @export
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
RprimerProfile <- function(x, ...) {
    methods::new("RprimerProfile", SummarizedExperiment(list(x = x), ...))
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
    # if (any(
    # grepl(paste0("[^", dnaBases, "]"),
    # row.names(SummarizedExperiment::assay(object))))
    # ) {
    #  msg <- c(msg, "The object has invalid rownames.")
    # }
    if (is.null(msg)) {
        TRUE
    } else {
        msg
    }
})
