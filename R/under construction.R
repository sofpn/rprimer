
# tm match primer


# 3end comp


# cons prof mask + plot
## tester
## assay


#' @noRd
#'
#' @keywords internal
#'
#' @examples
#' data("exampleRprimerOligo")
#' data("exampleRprimerAlignment")
#' target <- exampleRprimerAlignment
#' x <- exampleRprimerOligo$sequence[[1]]
#' .identifyBindingRegion(x, target)
.identifyBindingRegion <- function(x, target) {
    x <- Biostrings::DNAStringSet(x)[[1]]
    target <- Biostrings::DNAStringSet(target)
    match <- Biostrings::vmatchPattern(x, target, max.mismatch = 5)
    match <- as.data.frame(match)
    start <- sort(table(match$start), decreasing = TRUE)[1]
    end <- sort(table(match$end), decreasing = TRUE)[1]
    c("start" = as.numeric(names(start)), "end" = as.numeric(names(end)))
}
