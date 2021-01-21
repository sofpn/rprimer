
# mixed oligos

#' Generate oligos of a specific length
#'
#' \code{.generateOligos()} is the first step of the oligo-design process.
#' It finds all possible oligos of a specific length from an
#' \code{RprimerProfile} object, and returns a list containing
#' start and end position,
#' length, IUPAC sequence (the oligo DNA sequence with wobble bases), degeneracy
#' (number of variants of each oligo), maximum gap frequency,
#' mean overall identity, and minimum 3'-end identity at
#' both forward and reverse direction.
#'
#' @param x An \code{RprimerProfile} object.
#'
#' @param lengthOligo
#'
#' @return A list.
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerProfile")
#' .generateOligos(exampleRprimerProfile, lengthOligo = 18)
.generateOligos <- function(x, lengthOligo = 20) {
    oligos <- list()
    oligos$majoritySequence <- .nmers(x$majority, lengthOligo)
    oligos$iupacSequence <- .nmers(x$iupac, lengthOligo)
    oligos$start <- seq_len(nrow(oligos$iupacSequence)) + min(x$position) - 1
    oligos$end <- seq_len(
        nrow(oligos$iupacSequence)
    ) + lengthOligo - 1 + min(x$position) - 1
    oligos$length <- rep(lengthOligo, nrow(oligos$iupacSequence))
    oligos$degeneracy <- apply(oligos$iupacSequence, 1, .countDegeneracy)
    oligos$gapFrequency <- apply(.nmers(x$gaps, lengthOligo), 1, max)
    oligos$identity <- .nmers(x$identity, lengthOligo)
    oligos$coverage <- .nmers(x$coverage, lengthOligo)
    oligos$endCoverageFwd <- apply(
        oligos$coverage[
            , (ncol(oligos$coverage) - 6):ncol(oligos$coverage)],
        1, min
    )
    oligos$endCoverageRev <- apply(oligos$coverage[, seq_len(6)], 1, min)
    oligos$coverage <- rowMeans(oligos$coverage)
    oligos$roiStart <- rep(
        min(x$position, na.rm = TRUE), nrow(oligos$iupacSequence)
    )
    oligos$roiEnd <- rep(
        max(x$position, na.rm = TRUE), nrow(oligos$iupacSequence)
    )
    oligos
}



#' Combine majority and IUPAC
#' split in halves
#'
.splitAndPaste <- function(first, second) {
    first <- first[, 1:(ncol(first) / 2), drop = FALSE]
    second <- second[, (ncol(second) / 2 + 1):ncol(second), drop = FALSE] ########################################################
    cbind(first, second)
}

.addMixedPrimers <- function(x) {
    x$mixedFwd <- .splitAndPaste(x$majoritySequence, x$iupacSequence)
    x$mixedRev <- .splitAndPaste(x$iupacSequence, x$majoritySequence)
    x$mixedRev <- .reverseComplement(x$mixedRev)
    x$degeneracyMixedFwd <- apply(x$mixedFwd, 1, .countDegeneracy)
    x$degeneracyMixedRev <- apply(x$mixedRev, 1, .countDegeneracy)
    x

    # half identity and half coverage

    # coverage will also be different
}
