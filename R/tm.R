#' Split a DNA sequence into nearest neighbors
#'
#' @param x A vector with a DNA sequence.
#'
#' @return A vector with the nearest neighbors of \code{x}.
#'
#' @keywords internal
#'
#' @noRd
.nn <- function(x) {
    from <- (seq_along(x) - 1)[-1]
    to <- seq_along(x)[-1]
    vapply(seq_along(from), function(i) {
        paste(x[from[i]:to[i]], collapse = "")
    }, character(1))
}

#' Calculate stack values for dH or dS of nearest neighbors
#'
#' @param x A matrix with nearest-neighbor pairs
#' of DNA sequences.
#'
#' @param table The lookup table that should be used.
#' Either 'dH' (entropy) or 'dS' (enthalpy).
#'
#' @return The corresponding values for dH or dS (in cal/M). A numeric matrix.
#'
#' @keywords internal
#'
#' @noRd
.stack <- function(x, table = "dH") {
    selected <- lookup$nn[[table]]
    stack <- matrix(nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
    stack[] <- selected[match(x, lookup$nn$bases)]
    stack
}

#' Calculate initiation values for dH or dS of nearest neighbors
#'
#' @param x A matrix with nearest-neighbor pairs
#' of DNA sequences.
#'
#' @param table The lookup table that should be used.
#' Either 'dH' (entropy) or 'dS' (enthalpy).
#'
#' @return The corresponding values for dH or dS (in cal/M). A numeric vector.
#'
#' @keywords internal
#'
#' @noRd
.initiate <- function(x, table = "dH") {
    initiation <- lookup$nn[lookup$nn$bases == "Initiation", ][[table]]
    penalty <- lookup$nn[lookup$nn$bases == "AT_penalty", ][[table]]
    penaltyFirst <- ifelse(x[, 1] == "AT" | x[, 1] == "TA", penalty, 0)
    penaltyLast <- ifelse(
        x[, ncol(x)] == "AT" | x[, ncol(x)] == "TA", penalty, 0
    )
    initiation + penaltyFirst + penaltyLast
}

#' Get delta H and delta S values for Tm calculation
#'
#' \code{.tmParameters()} calculates and sums delta H and delta S for one or
#' more perfectly matching DNA duplexes (oligo-target duplexes),
#' using the nearest neighbor method. It also returns the number
#' of phosphates and the sodium ion concentration. For Tm calculation.
#'
#' @param x A matrix with DNA sequences.
#'
#' @param concNa Sodium ion concentration in M.
#'
#' @return A numeric matrix.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' .tmParameters(matrix(c("A", "G", "T", "T", "C", "G", "G", "T", "C", "G")))
.tmParameters <- function(x, concNa = 0.05) {
    nn <- t(apply(x, 1, .nn))
    dhStack <- rowSums(.stack(nn, "dH"))
    dsStack <- rowSums(.stack(nn, "dS"))
    dhInit <- .initiate(nn, "dH")
    dsInit <- .initiate(nn, "dS")
    sumdH <- dhStack + dhInit
    sumdS <- dsStack + dsInit
    n <- dim(x)[2] - 1 ## Number of phosphates
    m <- matrix(
        c(sumdH, sumdS, rep(n, nrow(nn)), rep(concNa, nrow(nn))),
        ncol = 4
    )
    rownames(m) <- rownames(nn)
    colnames(m) <- c("sumdH", "sumdS", "n", "concNa")
    m
}

#' Calculate melting temperature
#'
#' @param x An output from \code{.tmParameters()}
#'
#' @param concOligo Oligo concentration in nM.
#'
#' @keywords internal
#'
#' @noRd
.tm <- function(x, concOligo = 250) {
    concOligo <- concOligo * 10^(-9)
    tm <- x["sumdH"] / (x["sumdS"] + 0.368 * x["n"] * log(x["concNa"]) + 1.987 * log(concOligo)) - 273.15
    unname(tm)
}
