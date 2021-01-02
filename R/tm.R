#' Calculate melting temperature
#'
#' \code{.tm()} calculates the melting temperature of one or
#' more perfectly matching DNA duplexes (i.e. oligo-target duplexes),
#' using the nearest neighbor method.
#'
#' @param x
#' A character vector or matrix, where each row corresponds to a specific oligo.
#'
#' @param concOligo
#' Oligo concentration in nM. A number
#' [20, 2000] Defaults to 250.
#'
#' @param concNa
#' Sodium ion concentration in the PCR reaction in M.
#' A number [0.01, 1]. Defaults to 0.05 M (50 mM).
#'
#' @return The melting temperature(s) of \code{x}. A numeric vector.
#'
#' @section Details:
#'
#' Melting temperatures are calculated using SantaLucia's nearest-neighbor
#' method, with the following assumptions:
#'
#' \itemize{
#'   \item Oligos are not expected to be self-complementary, and hence no
#'   symmetry correction is done.
#'   \item The oligo concentration is assumed to be much higher
#'   than the target concentration.
#' }
#'
#' @examples
#' .tm(c("A", "G", "T", "T", "C", "G", "G", "T", "C", "G"))
#'
#' @references
#' SantaLucia Jr, J., & Hicks, D. (2004).
#' The thermodynamics of DNA structural motifs.
#' Annu. Rev. Biophys. Biomol. Struct., 33, 415-440.
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
    if (!is.matrix(x)) x <- t(matrix(x))
    nn <- t(apply(x, 1, .nnSplit))
    dhStack <- rowSums(.getStackValue(nn, "dH"))
    dsStack <- rowSums(.getStackValue(nn, "dS"))
    dhInit <- .getInitiationValue(nn, "dH")
    dsInit <- .getInitiationValue(nn, "dS")
    sumdH <- dhStack + dhInit
    sumdS <- dsStack + dsInit
    N <- dim(x)[2] - 1 # Number of phosphates
    sumdH / (sumdS + 0.368 * N * log(concNa) + 1.987 * log(concOligo)) - 273.15
}

# Helpers ======================================================================

#' Split a DNA sequence into nearest neighbors
#'
#' \code{.nnSplit()} splits a DNA sequence into nearest neighbors
#' (for calculation of deltaS, deltaH, deltaG and Tm).
#' Helper function to \code{.tm()}.
#'
#' @param x a matrix with DNA sequences.
#'
#' @return The nearest neighbors of x (a matrix).
#'
#' @keywords internal
#'
#' @noRd
.nnSplit <- function(x) {
    from <- (seq_along(x) - 1)[-1]
    to <- seq_along(x)[-1]
    vapply(seq_along(from), function(i) {
        paste(x[from[i]:to[i]], collapse = "")
    }, character(1))
}

#' Calculate stack values for dH or dS of nearest neighbors
#'
#' \code{.getStackValue()} finds the corresponding dH or dS values
#' for nearest-neighbor pairs. Helper function to \code{.tm()}.
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
    tableValues <- lookup$nn[[table]]
    stack <- tableValues[match(x, lookup$nn$bases)]
    if (!is.null(ncol(x))) {
        stack <- matrix(stack, ncol = ncol(x), byrow = FALSE)
    }
    stack
}

#' Calculate initiation values for dH or dS of nearest neighbors
#'
#' \code{.getInitiationValue()} finds the corresponding dH or dS initiation
#' values. Helper function to \code{.tm()}.
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
    initiation <- lookup$nn[lookup$nn$bases == "Initiation", ][[table]]
    penalty <- lookup$nn[lookup$nn$bases == "AT_penalty", ][[table]]
    penaltyFirst <- ifelse(x[, 1] == "AT" | x[, 1] == "TA", penalty, 0)
    penaltyLast <- ifelse(
        x[, ncol(x)] == "AT" | x[, ncol(x)] == "TA", penalty, 0
    )
    initiation + penaltyFirst + penaltyLast
}

# .adjustOligoConc <- function(x, oldOligoConc, newOligoConc) {
#}
#  tm =   sumdH / (sumdS + 0.368 * N * log(concNa) + 1.987 * log(500)) - 273.15



