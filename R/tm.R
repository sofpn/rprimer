#' Split a DNA sequence into nearest neighbors
#'
#' @param x A vector with a DNA sequence.
#'
#' @return A vector with the nearest neighbors of \code{x}.
#'
#' @examples
#' .nn(c("A", "G", "C", "T"))
.nn <- function(x) {
    from <- (seq_along(x) - 1)[-1]
    to <- seq_along(x)[-1]
    vapply(seq_along(from), function(i) {
        paste(x[from[[i]]:to[[i]]], collapse = "")
    }, character(1L))
}

#' Calculate stack values for dH or dS of nearest neighbors
#'
#' @examples
#' x <- t(matrix(.nn(c("A", "C", "G"))))
#' .stack(x)
.stack <- function(x, table = "dH") {
    selected <- lookup$nn[[table]]
    stack <- matrix(nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
    stack[] <- selected[match(x, lookup$nn$bases)]
    stack
}

#' Calculate initiation values for dH or dS of nearest neighbors
#'
#' @examples
#' x <- t(matrix(.nn(c("A", "C", "G"))))
#' .initiate(x)
.initiate <- function(x, table = "dH") {
    initiation <- lookup$nn[lookup$nn$bases == "Initiation", ][[table]]
    penalty <- lookup$nn[lookup$nn$bases == "AT_penalty", ][[table]]
    penaltyFirst <- ifelse(grepl("^(A|T)", x[, 1]), penalty, 0)
    penaltyLast <- ifelse(grepl("(A|T)$", x[, ncol(x)]), penalty, 0)
    initiation + penaltyFirst + penaltyLast
}

#' @examples
#' .tmParameters(t(matrix(c("A", "G", "T", "T", "C", "G", "G", "T", "C", "G"))))
.tmParameters <- function(x, concNa = 0.05) {
    nn <- t(apply(x, 1, .nn))
    dhStack <- rowSums(.stack(nn, "dH"))
    dsStack <- rowSums(.stack(nn, "dS"))
    dhInit <- .initiate(nn, "dH")
    dsInit <- .initiate(nn, "dS")
    sumdH <- dhStack + dhInit
    sumdS <- dsStack + dsInit
    n <- dim(x)[[2]] - 1 ## Number of phosphates
    sumdS <- sumdS + 0.368 * n * log(concNa) ## Salt correction of deltaS
    m <- matrix(c(sumdH, sumdS), ncol = 2)
    rownames(m) <- rownames(nn)
    colnames(m) <- c("sumdH", "sumdS")
    m
}

#' @examples
#' x <- .tmParameters(t(matrix(c("A", "G", "T", "T", "C", "G", "G", "T", "C"))))
#' .tm(x)
.tm <- function(x, concOligo = 250) {
    concOligo <- concOligo * 1e-9
    tm <- (x[, "sumdH"] * 1000) / (x[, "sumdS"] + 1.987 * log(concOligo)) - 273.15
    names(tm) <- rownames(x)
    tm
}
