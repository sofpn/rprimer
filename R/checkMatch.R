#' Check how oligos and assays match to target sequences
#'
#' \code{checkMatch()} checks how well the oligos within an \code{RprimerOligo}
#' or \code{RprimerAssay} object match with their intended target sequences.
#' The output gives information on both the proportion and names of target
#' sequences that match perfectly as well as with one, two, three or four or
#' more mismatches to the oligo within the intended oligo binding region
#' (i.e., on target match).
#' It also reports the proportion of target sequences that matches to the
#' oligo with no more than four mismatches within all other regions in
#' the alignment (off target match).
#' Ambiguous bases and gaps in the
#' target sequences will be identified as mismatches.
#' The function is a wrapper to \code{Biostrings::vcountPDict()}
#' (Pages et al., 2020).
#'
#' Note that the output does not say anything about the type,
#' position or severity of the mismatches.
#'
#' @param x
#' An \code{RprimerOligo} or \code{RprimerAssay} object.
#'
#' @param target
#' A \code{Biostrings::DNAMultipleAlignment} alignment with
#' intended target sequences. Note that it must be
#' same alignment that was used for generating the oligos/assays in \code{x}.
#'
#' @return
#' An \code{RprimerMatchOligo} or \code{RprimerMatchAssay} object, depending
#' on whether an \code{RprimerOligo} or \code{RprimerAssay} object was used
#' as input.
#'
#' \code{RprimerMatchOligo} objects contain the following information:
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
#'  \code{RprimerMatchAssat} objects contain the following information:
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
#' If the input assay contains probes, the following information
#' is also added:
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
#' @references
#' Pages, H., Aboyoun, P., Gentleman R., and DebRoy S. (2020). Biostrings:
#' Efficient manipulation of biological strings. R package version
#' 2.57.2.
setGeneric("checkMatch", function(x, target) standardGeneric("checkMatch"))

# Methods ======================================================================

#' @describeIn checkMatch
#'
#' @export
#' @examples
#' #### RprimerOligo objects
#'
#' data("exampleRprimerOligo")
#' data("exampleRprimerAlignment")
#'
#' x <- exampleRprimerOligo[1:2, ]
#' target <- exampleRprimerAlignment
#'
#' checkMatch(x, target)
#'
setMethod("checkMatch", "RprimerOligo", function(x, target) {
    if (!methods::is(target, "DNAMultipleAlignment")) {
        stop("'target' must be a DNAMultipleAlignment object.", call. = FALSE)
    }
    match <- .checkMatchOligo(x, target)
    RprimerMatchOligo(match)
})

#' Check match of an RprimerAssay object (method)
#'
#' @describeIn checkMatch
#'
#' @export
#'
#' @examples
#' #### RprimerAssay objects
#'
#' data("exampleRprimerAssay")
#' data("exampleRprimerAlignment")
#'
#' x <- exampleRprimerAssay[1:2, ]
#' target <- exampleRprimerAlignment
#'
#' checkMatch(x, target)
#'
setMethod("checkMatch", "RprimerAssay", function(x, target) {
    if (!methods::is(target, "DNAMultipleAlignment")) {
        stop("'target' must be a DNAMultipleAlignment object.", call. = FALSE)
    }
    match <- .checkMatchAssay(x, target)
    RprimerMatchAssay(match)
})

# Helpers ======================================================================

#' @noRd
#'
#' @examples
#' data("exampleRprimerAlignment")
#' .maskRange(from = 1, to = 10, exampleRprimerAlignment, invert = TRUE)
.maskRange <- function(from, to, target, ...) {
    selection <- target
    Biostrings::colmask(selection, ...) <- IRanges::IRanges(
        start = from, end = to
    )
    Biostrings::DNAStringSet(selection)
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerAlignment")
#' data("exampleRprimerOligo")
#' x <- exampleRprimerOligo$sequence[[1]]
#' x <- Biostrings::DNAStringSet(x)
#' target <- Biostrings::DNAStringSet(exampleRprimerAlignment)
#' .getMatchIndex(x, target)
.getMatchIndex <- function(x, target) {
    res <- lapply(seq(0, 3), function(i) {
        result <- Biostrings::vcountPDict(
            x, target,
            max.mismatch = i
        )
        which(colSums(result) == 0)
    })
    res[-length(res)] <- lapply(seq(1, length(res) - 1), function(i) {
        setdiff(res[[i]], res[[i + 1]])
    })
    res[[length(res) + 1]] <- setdiff(seq_along(target), unlist(res))
    res[c(length(res), seq_len(length(res) - 1))]
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerAlignment")
#' target <- Biostrings::DNAStringSet(exampleRprimerAlignment)
#' .getSequenceNames(1:2, target)
.getSequenceNames <- function(x, target) {
    names(target) <- sub(" .*", "", names(target))
    lapply(x, function(i) names(target)[i])
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerAlignment")
#' data("exampleRprimerOligo")
#' target <- Biostrings::DNAStringSet(exampleRprimerAlignment)
#' x <- exampleRprimerOligo$sequence[[1]]
#' x <- Biostrings::DNAStringSet(x)
#' .getMatchIndexOffTarget(x, target)
.getMatchIndexOffTarget <- function(x, target, max.mismatch = 4) {
    result <- Biostrings::vcountPDict(x, target, max.mismatch = max.mismatch)
    which(colSums(result) > 0)
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerAlignment")
#' data("exampleRprimerOligo")
#' target <- Biostrings::DNAStringSet(exampleRprimerAlignment)
#' x <- exampleRprimerOligo$sequence[[1]]
#' x <- Biostrings::DNAStringSet(x)
#' .getMatchProportion(x, target)
.getMatchProportion <- function(x, target) {
    matching <- .getMatchIndex(x, target)
    sequenceNames <- .getSequenceNames(matching, target)
    matching <- lapply(matching, function(x) length(x) / length(target))
    matching <- do.call("cbind.data.frame", matching)
    sequenceNames <- as.data.frame(t(do.call("cbind", list(sequenceNames))))
    matching <- cbind(matching, sequenceNames)
    matching <- matching[c(1, 6, 2, 7, 3, 8, 4, 9, 5, 10)]
    names(matching) <- c(
        "perfectMatch", "idPerfectMatch",
        "oneMismatch", "idOneMismatch",
        "twoMismatches", "idTwoMismatches",
        "threeMismatches", "idThreeMismatches",
        "fourOrMoreMismatches", "idFourOrMoreMmismatches"
    )
    matching
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerAlignment")
#' data("exampleRprimerOligo")
#' target <- Biostrings::DNAStringSet(exampleRprimerAlignment)
#' x <- exampleRprimerOligo$sequence[[1]]
#' x <- Biostrings::DNAStringSet(x)
#' .getMatchProportionOffTarget(x, target)
.getMatchProportionOffTarget <- function(x, target) {
    matching <- .getMatchIndexOffTarget(x, target)
    sequenceNames <- .getSequenceNames(list(matching), target)
    matching <- length(matching) / length(target)
    matching <- do.call("cbind.data.frame", list(matching))
    sequenceNames <- as.data.frame(t(do.call("cbind", list(sequenceNames))))
    matching <- cbind(matching, sequenceNames)
    names(matching) <- c("offTargetMatch", "idOffTargetMatch")
    matching
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerAlignment")
#' data("exampleRprimerOligo")
#' x <- exampleRprimerOligo[1:2, ]
#' target <- exampleRprimerAlignment
#' .checkMatchOligo(x, target)
.checkMatchOligo <- function(x, target) {
    onTarget <- lapply(seq_len(nrow(x)), function(i) {
        target <- .maskRange(
            x$start[[i]], x$end[[i]], target,
            invert = TRUE
        )
        check <- Biostrings::DNAStringSet(x$sequence[[i]])
        .getMatchProportion(check, target)
    })
    onTarget <- do.call("rbind", onTarget)
    offTarget <- lapply(seq_len(nrow(x)), function(i) {
        target <- .maskRange(
            x$start[[i]], x$end[[i]], target,
            invert = FALSE
        )
        check <- Biostrings::DNAStringSet(x$sequence[[i]])
        .getMatchProportionOffTarget(check, target)
    })
    offTarget <- do.call("rbind", offTarget)
    cbind("iupacSequence" = x$iupacSequence, onTarget, offTarget)
}

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
    fwd <- .checkMatchOligo(oligos$fwd, target)
    names(fwd) <- paste0(names(fwd), "Fwd")
    rev <- .checkMatchOligo(oligos$rev, target)
    names(rev) <- paste0(names(rev), "Rev")
    all <- cbind(fwd, rev)
    if (any(grepl("Pr$", names(x)))) {
        pr <- .checkMatchOligo(oligos$pr, target)
        names(pr) <- paste0(names(pr), "Pr")
        all <- cbind(all, pr)
    }
    all
}
