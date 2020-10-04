
#  New class ==================================================================

#' @export
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setClass("RprimerProperties", contains = "SummarizedExperiment")

# Constructor =================================================================

#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
RprimerProperties <- function(x, ...) {
  methods::new("RprimerProperties", SummarizedExperiment(list(x = x), ...))
}

# Validity method =============================================================

S4Vectors::setValidity2("RprimerProperties", function(object) {
  msg <- NULL

  if (is.null(msg)) {
    TRUE
  } else {
    msg
  }
})
