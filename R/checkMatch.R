#' Check how oligos match to their target sequences (generic)
#'
#' @param x
#' An \code{RprimerOligo} or \code{RprimerAssay} object.
#'
#' @param target
#' The \code{Biostrings::DNAMultipleAlignment} object used for oligo/assay
#' design.
#'
#' ## Write details and limitations
#' - only taerg region selected, seq w amb bases rmoved, no pos
#'
#'  the check will performed at all regions within the alignment except
#'  for the intended binding region (to detect any potential issues with off
#'  target binding).
#'  performed at the intended binding region (all other regions will be
#'  excluded to not generate false ). Here, to avoid false negative results,
#'  sequences with ambiguous bases at the oligo binding region are removed.
#'
#' The function is a wrapper to \code{Biostrings::vcountPDict()}
#' (Pages et al., 2020)
#'
#' @return
#' An \code{RprimerMatchOligo} or \code{RprimerMatchAssay} object. See details
#' below.
#'
#' @export
#'
#' @references
#' Pages, H., Aboyoun, P., Gentleman R., and DebRoy S. (2020). Biostrings:
#' Efficient manipulation of biological strings. R package version
#' 2.57.2.
setGeneric("checkMatch", function(x, target) standardGeneric("checkMatch"))

# Methods ======================================================================

#' Check match of an RprimerOligo object (method)
#'
#' @describeIn checkMatch
#'
#' @export
#'
#' @section Output (for oligos):
#'
#' The output contains the following information:
#'
#' \describe{
#'   \item{iupacSequence}{The oligo sequence in IUPAC format.}
#'   \item{perfectMatch}{Proportion of target sequences that matches perfectly
#'   to the oligo within the intended binding region.}
#'   \item{idPerfectMatch}{Names of all sequences that matches perfectly.}
#'   \item{oneMismatch}{Proportion of target sequences with one mismatch to
#'   the oligo
#'   within the intended binding region.}
#'   \item{idOneMismatch}{Names of all sequences that matches with one
#'   mismatch.}
#'   \item{twoMismatches}{Proportion of target sequences with two mismatches
#'   to the oligo
#'   within the intended binding region.}
#'   \item{idTwoMismatches}{Names of all sequences that matches with two
#'   mismatches.}
#'   \item{threeMismatches}{Proportion of target sequences with
#'   three mismatches to the oligo
#'   within the intended binding region.}
#'   \item{idThreeMismatches}{Names of all sequences that matches with three
#'   mismatches.}
#'   \item{fourOrMoreMismatches}{Proportion of target sequences with
#'   four or more mismatches to the oligo
#'   within the intended binding region.}
#'   \item{idFourOrMoreMismatches}{Names of all sequences that matches with four
#'   or more mismatches.}
#'   \item{offTargetMatch}{Proportion of target sequences that matches
#'   to the oligo with
#'   no more than four mismatches to all other regions within the alignment.}
#'   \item{idOffTargetMatch}{Names of all sequences that matches with
#'   no more than four mismatches to all other regions within the alignment.}
#'  }
#'
#'
#' @examples
#' data("exampleRprimerOligo")
#' data("exampleRprimerAlignment")
#' x <- exampleRprimerOligo[1:2, ]
#' target <- exampleRprimerAlignment
#' checkMatch(x, target)
setMethod("checkMatch", "RprimerOligo", function(x, target) {
    if (!methods::is(target, "DNAMultipleAlignment")) {
        stop("'target' must be a DNAMultipleAlignment object.", call. = FALSE)
    }
    match <- .checkMatchOligos(x, target)
    RprimerMatchOligo(match)
})

#' Check match of an RprimerAssay object (method)
#'
#' @describeIn checkMatch
#'
#' @section Output (for assays):
#'
#' The output contains the following information:
#'
#' \describe{
#'   \item{iupacSequenceFwd}{The forward primer sequence in IUPAC format.}
#'   \item{perfectMatchFwd}{Proportion of target sequences that matches
#'   perfectly
#'   with the forward primer withing the intended binding region.}
#'   \item{idPerfectMatchFwd}{Names of all sequences that matches perfectly.}
#'   \item{oneMismatchFwd}{Proportion of target sequences with one mismatch to
#'   the forward primer
#'   within the intended binding region.}
#'   \item{idOneMismatchFwd}{Names of all sequences that matches with one
#'   mismatch.}
#'   \item{twoMismatchesFwd}{Proportion of target sequences with two mismatches
#'   to the forward primer
#'   within the intended binding region.}
#'   \item{idTwoMismatchesFwd}{Names of all sequences that matches with two
#'   mismatches.}
#'   \item{threeMismatchesFwd}{Proportion of target sequences with
#'   three mismatches to the forward primer
#'   within the intended binding region.}
#'   \item{idThreeMismatchesFwd}{Names of all sequences that matches with three
#'   mismatches.}
#'   \item{fourOrMoreMismatchesFwd}{Proportion of target sequences with
#'   four or more mismatches to the forward primer
#'   within the intended binding region.}
#'   \item{idFourOrMoreMismatchesFwd}{Names of all sequences that matches with
#'   four or more mismatches.}
#'   \item{offTargetMatchFwd}{Proportion of target sequences that matches
#'   to the forward primer with
#'   no more than four mismatches to all other regions within
#'   the target alignment.}
#'   \item{idOffTargetMatchFwd}{Names of all target sequences that matches
#'   to the forward primer with
#'   no more than four mismatches to all other regions within
#'   the target alignment.}
#'   \item{iupacSequenceRev}{The reverse primer sequence in IUPAC format.}
#'   \item{perfectMatchRev}{Proportion of target sequences that matches
#'   perfectly
#'   with the reverse primer withing the intended binding region.}
#'   \item{idPerfectMatchRev}{Names of all sequences that matches perfectly.}
#'   \item{oneMismatchRev}{Proportion of target sequences with one mismatch to
#'   the reverse primer
#'   within the intended binding region.}
#'   \item{idOneMismatchRev}{Names of all sequences that matches with one
#'   mismatch.}
#'   \item{twoMismatchesRev}{Proportion of target sequences with two mismatches
#'   to the reverse primer
#'   within the intended binding region.}
#'   \item{idTwoMismatchesRev}{Names of all sequences that matches with two
#'   mismatches.}
#'   \item{threeMismatchesRev}{Proportion of target sequences with
#'   three mismatches to the reverse primer
#'   within the intended binding region.}
#'   \item{idThreeMismatchesRev}{Names of all sequences that matches with three
#'   mismatches.}
#'   \item{fourOrMoreMismatchesRev}{Proportion of target sequences with
#'   four or more mismatches to the reverse primer
#'   within the intended binding region.}
#'   \item{idFourOrMoreMismatchesRev}{Names of all sequences that matches with
#'   four or more mismatches.}
#'   \item{offTargetMatchRev}{Proportion of target sequences that matches
#'   to the reverse primer with
#'   no more than four mismatches to all other regions within
#'   the target alignment.}
#'   \item{idOffTargetMatchRev}{target sequences that matches
#'   to the reverse primer with
#'   no more than four mismatches to all other regions within
#'   the target alignment.}
#'   \item{iupacSequencePr}{The probe sequence in IUPAC format.}
#'
#'   If the input assay contains probes, the following information
#'   is also added:
#'
#'   \item{perfectMatchPr}{Proportion of target sequences that matches perfectly
#'   with the probe withing the intended binding region.}
#'   \item{idPerfectMatchPr}{Names of all sequences that matches perfectly.}
#'   \item{oneMismatchPr}{Proportion of target sequences with one mismatch to
#'   the probe
#'   within the intended binding region.}
#'   \item{idOneMismatchPr}{Names of all sequences that matches with one
#'   mismatch.}
#'   \item{twoMismatchesPr}{Proportion of target sequences with two mismatches
#'   to the probe
#'   within the intended binding region.}
#'   \item{idTwoMismatchesPr}{Names of all sequences that matches with two
#'   mismatches.}
#'   \item{threeMismatchesPr}{Proportion of target sequences with
#'   three mismatches to the probe
#'   within the intended binding region.}
#'   \item{idThreeMismatchesPr}{Names of all sequences that matches with three
#'   mismatches.}
#'   \item{fourOrMoreMismatchesPr}{Proportion of target sequences with
#'   four or more mismatches to the probe
#'   within the intended binding region.}
#'   \item{idFourOrMoreMismatchesPr}{Names of all sequences that matches with
#'   four or more mismatches.}
#'   \item{offTargetMatchPr}{Proportion of target sequences that matches
#'   to the probe with
#'   no more than four mismatches to all other regions within
#'   the target alignment.}
#'   \item{idOffTargetMatchPr}{Names of all target sequences that matches
#'   to the probe with
#'   no more than four mismatches to all other regions within
#'   the target alignment.}
#'  }
#'
#' @export
#'
#' @examples
#' data("exampleRprimerAssay")
#' data("exampleRprimerAlignment")
#' x <- exampleRprimerAssay[1:2, ]
#' target <- exampleRprimerAlignment
#' checkMatch(x, target)
setMethod("checkMatch", "RprimerAssay", function(x, target) {
    if (!methods::is(target, "DNAMultipleAlignment")) {
        stop("'target' must be a DNAMultipleAlignment object.", call. = FALSE)
    }
    match <- .checkMatchAssay(x, target)
    RprimerMatchAssay(match)
})

# Helpers ======================================================================

.extractRange <- function(from, to, target, excludeAmbiguous = TRUE, ...) {
    selection <- target
    Biostrings::colmask(selection, ...) <- IRanges::IRanges(
        start = from, end = to
    )
    selection <- Biostrings::DNAStringSet(selection)
    if (excludeAmbiguous) selection <- ShortRead::clean(selection)
    selection
}

.getMatchIndex <- function(x, target) {
    res <- lapply(seq(0, 3), function(i) {
        result <- Biostrings::vcountPDict(
            x, target, max.mismatch = i
        )
        which(colSums(result) == 0)
    })
    res[-length(res)] <- lapply(seq(1, length(res) - 1), function(i) {
        setdiff(res[[i]], res[[i + 1]])
    })
    res[[length(res) + 1]] <- setdiff(seq_along(target), unlist(res))
    res[c(length(res), seq_len(length(res) - 1))]
}

.getSequenceNames <- function(x, target) {
    names(target) <- sub(" .*", "", names(target))
    lapply(x, function(i) names(target)[i])
}

.getMatchIndexOffTarget <- function(x, target) {
    result <- Biostrings::vcountPDict(x, target, max.mismatch = 4)
    which(colSums(result) > 0)
}

.getMatchProportion <- function(x, target) {
    matching <- .getMatchIndex(x, target)
    sequenceNames <- .getSequenceNames(matching, target)
    matching <- lapply(matching, function(x) length (x) / length(target))
    matching <- do.call("cbind.data.frame", matching)
    sequenceNames <- as.data.frame(t(do.call("cbind", list(sequenceNames))))
    matching <- cbind(matching, sequenceNames)
    matching <- matching[c(1, 6, 2, 7, 3, 8, 4, 9, 5, 10)]
    names(matching) <- c(
        "perfectMatch", "idPerfectMatch",
        "oneMismatch", "idOneMismatch",
        "twoMismatches", "idTwoMmismatches",
        "threeMismatches", "idThreeMmismatches",
        "fourOrMoreMismatches", "idFourOrMoreMmismatches"
    )
    matching
}

.getMatchProportionOffTarget <- function(x, target) {
    matching <- .getMatchIndexOffTarget(x, target)
    sequenceNames <- .getSequenceNames(list(matching), target)
    matching <- length(matching) /  length(target)
    matching <- do.call("cbind.data.frame", list(matching))
    sequenceNames <- as.data.frame(t(do.call("cbind", list(sequenceNames))))
    matching <- cbind(matching, sequenceNames)
    names(matching) <- c("offTargetMatch", "idOffTargetMatch")
    matching
}

.checkMatchOligos <- function(x, target) {
    onTarget <- lapply(seq_len(nrow(x)), function(i) {
        target <- .extractRange(
            x$start[[i]], x$end[[i]], target, invert = TRUE
        )
        check <- Biostrings::DNAStringSet(x$sequence[[i]])
        .getMatchProportion(check, target)
    })
    onTarget <- do.call("rbind", onTarget)
    offTarget <- lapply(seq_len(nrow(x)), function(i) {
        target <- .extractRange(
            x$start[[i]], x$end[[i]], target,
            excludeAmbiguous = FALSE, invert = FALSE
        )
        check <- Biostrings::DNAStringSet(x$sequence[[i]])
        .getMatchProportionOffTarget(check, target)
    })
    offTarget <- do.call("rbind", offTarget)
    cbind("iupacSequence" = x$iupacSequence, onTarget, offTarget)
}

#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerAssay")
#' x <- eampleRprimerAssay[1:4, ]
#' .convertAssay(x)
.convertAssay <- function(x) {
    x <- as.data.frame(x)
    fwd <- x[, grepl("Fwd", names(x))]
    names(fwd) <- gsub("Fwd", "", names(fwd))
    rev <- x[, grepl("Rev", names(x))]
    names(rev) <- gsub("Rev", "", names(rev))
    rev$sequence <- lapply(rev$sequence, function(x) {
        x <- Biostrings::DNAStringSet(x)
        x <- Biostrings::reverseComplement(x)
        as.character(x)
    })
    pr <- x[, grepl("Pr", names(x))]
    names(pr) <- gsub("Pr", "", names(pr))
    list("fwd" = fwd, "rev" = rev, "pr" = pr)
}

#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' data("exampleRprimerAlignment")
#' data("exampleRprimerAssay")
#' x <- exampleRprimerAssay[1:2, ]
#' target <- exampleRprimerAlignment
#' .checkMatchAssay(x, target)
.checkMatchAssay <- function(x, target) {
    oligos <- .convertAssay(x)
    fwd <- .checkMatchOligos(oligos$fwd, target)
    names(fwd) <- paste0(names(fwd), "Fwd")
    rev <- .checkMatchOligos(oligos$rev, target)
    names(rev) <- paste0(names(rev), "Rev")
    all <- cbind(fwd, rev)
    if (any(grepl("Pr$", names(x)))) {
        pr <- .checkMatchOligos(oligos$pr, target)
        names(pr) <- paste0(names(pr), "Pr")
        all <- cbind(all, pr)
    }
    all
}
