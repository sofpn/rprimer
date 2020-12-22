#' Split a DNA sequence into nearest neighbors
#'
#' \code{.nnSplit()} splits a DNA sequence into nearest neighbors
#' (for calculation of deltaG, deltaH and Tm)
#'
#' @param x a DNA sequence with at least two bases, e.g. 'CTTA'
#' (a character vector of length one).
#'
#' @return The nearest neighbors of x (a character vector).
#'
#' @keywords internal
#'
#' @noRd
.nnSplit <- function(x) {
    x <- .splitSequence(x)
    from <- (seq_along(x) - 1)[-1]
    to <- seq_along(x)[-1]
    nn <- purrr::map2_chr(from, to, function(i, j) {
        paste(x[i:j], collapse = "")
    })
    nn
}

#' Calculate stack values for dH or dS of nearest neighbors
#'
#' \code{.getStackValue()} finds the corresponding dH or dS values
#' for nearest-neighbor pairs.
#'
#' @param x A matrix or vector (of type character) with nearest-neighbor pairs
#' of DNA sequences (e.g. \code{c('CT', 'TT', 'TA')})
#'
#' @param table The lookup table that should be used.
#' Either 'dH' (entropy) or 'dS' (enthalpy).
#'
#' @return The corresponding values for dH or dS (in cal/M).
#'
#' @keywords internal
#'
#' @noRd
.getStackValue <- function(x, table = "dH") {
    tableValues <- rprimerGlobals$nnLookup[[table]]
    stack <- tableValues[match(x, rprimerGlobals$nnLookup$bases)]
    if (!is.null(ncol(x))) {
        stack <- matrix(stack, ncol = ncol(x), byrow = FALSE)
    }
    stack
}

#' Calculate initiation values for dH or dS of nearest neighbors
#'
#' \code{.getInitiationValue()} finds the corresponding dH or dS values
#' for nearest-neighbor pairs.
#'
#' @param x A matrix or vector (of type character) with nearest-neighbor pairs
#' of DNA sequences (e.g. \code{c('CT', 'TT', 'TA')})
#'
#' @param table The lookup table that should be used.
#' Either 'dH' (entropy) or 'dS' (enthalpy).
#'
#' @return The corresponding values for dH or dS (in cal/M).
#'
#' @keywords internal
#'
#' @noRd
.getInitiationValue <- function(x, table = "dH") {
    initValue <- rprimerGlobals$nnLookup[
        rprimerGlobals$nnLookup$bases == "Initiation",
    ][[table]]
    penaltyValue <- rprimerGlobals$nnLookup[
        rprimerGlobals$nnLookup$bases == "Terminal_AT_penalty",
    ][[table]]
    first <- x[, 1]
    last <- x[, ncol(x)]
    penaltyFirst <- ifelse(first == "AT" | first == "TA", penaltyValue, 0)
    penaltyLast <- ifelse(last == "AT" | last == "TA", penaltyValue, 0)
    value <- penaltyFirst + penaltyLast + initValue
    value
}

#' Melting temperature
#'
#' \code{.tm} calculates the melting temperature of one or
#' more perfectly matching DNA duplexes (i.e. oligo-target duplexes),
#' using the nearest neighbor method.
#'
#' @param oligos One or more DNA sequences (a character vector).
#'
#' @inheritParams getOligos
#'
#' @return The melting temperature(s) of x.
#'
#' @keywords internal
#'
#' @noRd
.tm <- function(x, concOligo = 500, concNa = 0.05) {
    if (concOligo < 20 || concOligo > 2000) {
        stop("'concOligo' must be from 20 nM to 2000 nM.", call. = FALSE)
    }
    if (concNa < 0.01 || concNa > 1) {
        stop("'concNa' must be from 0.01 to 1 M.", call. = FALSE)
    }
    concOligo <- concOligo * 10^(-9)
    x <- toupper(x)
    nn <- purrr::map(x, .nnSplit)
    oligoLength <- purrr::map_int(nn, length)
    if (length(unique(oligoLength)) != 1) {
        stop("All oligos must be of equal length.", call. = FALSE)
    }
    nn <- do.call("rbind", nn)
    dhStack <- rowSums(.getStackValue(nn, "dH"))
    dsStack <- rowSums(.getStackValue(nn, "dS"))
    dhInit <- .getInitiationValue(nn, "dH")
    dsInit <- .getInitiationValue(nn, "dS")
    sumdH <- dhStack + dhInit
    sumdS <- dsStack + dsInit
    N <- nchar(x[[1]]) - 1 # Number of phosphates
    tm <- sumdH / (sumdS + 0.368 * N * log(concNa) + 1.987 * log(concOligo))
    tm <- tm - 273.15
    tm
}
