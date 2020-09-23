
#  New class ==================================================================

#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.RprimerProfile <- setClass(
  "RprimerProfile",
  slots = representation(
    rowVec = "numeric",
    colVec = "numeric",
    rowToRowMat = "matrix",
    colToColMat = "matrix",
    rowToColMat = "matrix",
    colToRowMat = "matrix"
  ),
  contains = "SummarizedExperiment"
)

# Constrctor ==================================================================

#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
RprimerProfile <- function(
  rowVec = 0,
  colVec = 0,
  rowToRowMat = matrix(0,0,0),
  colToColMat = matrix(0,0,0),
  rowToColMat = matrix(0,0,0),
  colToRowMat = matrix(0,0,0),
  ...) {
  se <- SummarizedExperiment(...)
  .RprimerProfile(
    se,
    rowVec = rowVec,
    colVec = colVec,
    rowToRowMat = rowToRowMat,
    colToColMat = colToColMat,
    rowToColMat = rowToColMat,
    colToRowMat = colToRowMat
    )
}

# Getter methods for 1D data ==================================================

#' @export
setGeneric("rowVec", function(x, ...) standardGeneric("rowVec"))

#' @export
setGeneric("colVec", function(x, ...) standardGeneric("colVec"))

#' @export
setMethod("rowVec", "RprimerProfile", function(x, withDimnames=TRUE) {
  out <- x@rowVec
  if (withDimnames)
    names(out) <- rownames(x)
  out
})

#' @export
setMethod("colVec", "RprimerProfile", function(x, withDimnames=TRUE) {
  out <- x@colVec
  if (withDimnames)
    names(out) <- colnames(x)
  out
})

# Getter methods for 2D data ==================================================

#' @export
setGeneric("rowToRowMat", function(x, ...) standardGeneric("rowToRowMat"))

#' @export
setGeneric("colToColMat", function(x, ...) standardGeneric("colToColMat"))

#' @export
setGeneric("rowToColMat", function(x, ...) standardGeneric("rowToColMat"))

#' @export
setGeneric("colToRowMat", function(x, ...) standardGeneric("colToRowMat"))

#' @export
setMethod("rowToRowMat", "RprimerProfile", function(x, withDimnames=TRUE) {
  out <- x@rowToRowMat
  if (withDimnames)
    rownames(out) <- rownames(x)
  out
})

#' @export
setMethod("colToColMat", "RprimerProfile", function(x, withDimnames=TRUE) {
  out <- x@colToColMat
  if (withDimnames)
    colnames(out) <- colnames(x)
  out
})

#' @export
setMethod("rowToColMat", "RprimerProfile", function(x, withDimnames=TRUE) {
  out <- x@rowToColMat
  if (withDimnames)
    rownames(out) <- colnames(x)
  out
})

#' @export
setMethod("colToRowMat", "RprimerProfile", function(x, withDimnames=TRUE) {
  out <- x@colToRowMat
  if (withDimnames)
    colnames(out) <- rownames(x)
  out
})


# Validity method

#' @importFrom BiocGenerics NCOL NROW
setValidity2("RprimerProfile", function(object) {
  NR <- NROW(object)
  NC <- NCOL(object)
  msg <- NULL

  # 1D
  if (length(rowVec(object, withDimnames=FALSE)) != NR) {
    msg <- c(msg, "'rowVec' should have length equal to the number of rows")
  }
  if (length(colVec(object, withDimnames=FALSE)) != NC) {
    msg <- c(
      msg, "'colVec' should have length equal to the number of columns"
    )
  }

  # 2D
  if (NROW(rowToRowMat(object, withDimnames=FALSE)) != NR) {
    msg <- c(
      msg, "'nrow(rowToRowMat)' should be equal to the number of rows"
    )
  }
  if (NCOL(colToColMat(object, withDimnames=FALSE)) != NC) {
    msg <- c(
      msg, "'ncol(colToColMat)' should be equal to the number of columns"
    )
  }
  if (NROW(rowToColMat(object, withDimnames=FALSE)) != NC) {
    msg <- c(
      msg, "'nrow(rowToColMat)' should be equal to the number of columns"
    )
  }
  if (NCOL(colToRowMat(object, withDimnames=FALSE)) != NR) {
    msg <- c(
      msg, "'ncol(colToRowMat)' should be equal to the number of rows"
    )
  }

  if (length(msg)) {
    msg
  } else TRUE
})

# Show method =================================================================

#' @export
#' @importMethodsFrom SummarizedExperiment show
setMethod("show", "RprimerProfile", function(object) {
  callNextMethod()
  cat(
    "rowToRowMat has ", ncol(rowToRowMat(object)), " columns\n",
    "colToColMat has ", nrow(colToColMat(object)), " rows\n",
    "rowToColMat has ", ncol(rowToRowMat(object)), " columns\n",
    "colToRowMat has ", ncol(rowToRowMat(object)), " rows\n",
    sep=""
  )
})


























#S4Vectors::setValidity2("RprimerProfile", function(object) {
#  msg <- NULL
#  if (!is.matrix(SummarizedExperiment::assay(object))) {
#    msg <- c(msg, "The object is not a matrix.")
# }
#  if (!is.numeric(SummarizedExperiment::assay(object))) {
#    msg <- c(msg, "The object is not in numeric format.")
#  }
# if (min(SummarizedExperiment::assay(object)) < 0) {
#    msg <- c(msg, "The object contains values < 0.")
#  }
#  if (max(SummarizedExperiment::assay(object)) > 1) {
#    msg <- c(msg, "The object contains values > 1.")
#  }
#  if (rownames(SummarizedExperiment::assay(object)) = NULL) {
#    msg <- c(msg, "The object must have rownames.")
#  }
#  if (any(grepl(paste0("^", dnaBases), row.names(SummarizedExperiment::assay(object))))) {
#    msg <- c(msg, paste0(
#      "The object contains invalid rownames. \n Valid rownames are ",
#      dnaBases, "."
#    )
#  )
#  }
#  if (is.null(msg)) {
#    TRUE
#  } else {
#    msg
#  }
#})


# Best practice:

# AllClasses.R
# AllGenerics.R

