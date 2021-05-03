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
        paste(x[from[[i]]:to[[i]]], collapse = "")
    }, character(1L))
}

#' Calculate stack values for dH or dS of nearest neighbors
#'
#' @param x A matrix with nearest-neighbor pairs
#' of DNA sequences.
#'
#' @param table The lookup table that should be used.
#' Either 'dH' (entropy), 'dS' (enthalpy) or 'dG' (free energy)
#'
#' @return The corresponding values for dH, dS or dG (in cal/M).
#' A numeric matrix.
#'
#' @keywords internal
#'
#' @noRd
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
#' @param x A matrix with nearest-neighbor pairs
#' of DNA sequences.
#'
#' @param table The lookup table that should be used.
#' Either 'dH' (entropy), 'dS' (enthalpy) or 'dG' (free energy)
#'
#' @return The corresponding values for dH, dS or dG (in cal/M).
#' A numeric vector.
#'
#' @keywords internal
#'
#' @noRd
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

#' Calculate melting temperature
#'
#' @param x A named vector with tm parameters.
#'
#' @param concOligo Oligo concentration in nM.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' x <- .tmParameters(t(matrix(c("A", "G", "T", "T", "C", "G", "G", "T", "C"))))
#' .tm(x)
.tm <- function(x, concOligo = 250) {
    concOligo <- concOligo * 1e-9
    tm <- (x[, "sumdH"] * 1000) / (x[, "sumdS"] + 1.987 * log(concOligo)) - 273.15
    names(tm) <- rownames(x)
    tm
}
