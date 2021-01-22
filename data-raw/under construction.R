
#' Combine majority and IUPAC
#' split in halves
#'
.splitAndPaste <- function(first, second) {
    first <- first[, 1:(ncol(first) / 2), drop = FALSE]
    second <- second[, (ncol(second) / 2 + 1):ncol(second), drop = FALSE]
    cbind(first, second)
}

.addMixedPrimers <- function(x) {
    x$mixedFwd <- .splitAndPaste(x$majoritySequence, x$iupacSequence)
    x$mixedRev <- .splitAndPaste(x$iupacSequence, x$majoritySequence)
    x$degeneracyMixedFwd <- apply(x$mixedFwd, 1, .countDegeneracy)
    x$degeneracyMixedRev <- apply(x$mixedRev, 1, .countDegeneracy)
    x

    # half identity and half coverage
    # coverage will also be different
}
