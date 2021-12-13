# check only on unique seq varianrs


#' Check how oligos and assays match to their target sequences
#'
#' \code{checkMatch()} investigates how well oligos or assays match with their
#' intended target sequences
#' within a multiple DNA sequence alignment.
#'
#' The output provides information on the proportion and names of target
#' sequences
#' that match perfectly as well as with one, two, three, or four or more
#' mismatches to the oligo within the intended oligo binding region
#' in the input alignment (on-target match).
#' It also gives the proportion and names of target sequences that
#' match with a maximum of two mismatches to at least one sequence variant of
#' the oligo outside the
#' intended oligo binding region (off-target match).
#' The function is a wrapper to \code{Biostrings::vcountPDict()}
#' (Pages et al., 2020).
#'
#' @md
#'
#' @section Limitations:
#' There are a few limitations with this function, which is important to be
#' aware of:
#'
#' * False negatives or positives may occur due to
#' poorly aligned sequences
#' * The output does not tell which strand
#' (minus or plus) the oligo matches to.
#' This is important to consider when
#' assessing off-target matches to single-stranded
#' targets
#' * Ambiguous bases and gaps in the
#' target sequences are identified as mismatches
#' * The function checks strictly on- and off-target, and may therefore
#' miss off-target matches that partially overlap the intended target
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
#'   \item{offTargetMatch}{Proportion of target sequences with maximum two
#'   mismatches to at least one site outside the intended
#'   oligo binding region in the input alignment.
#'   }
#'   \item{idOffTargetMatch}{Names of all off-target matching sequences.}
#'
#' \code{RprimerMatchAssay} objects contain the following information:
#'
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
#'   \item{offTargetMatchFwd}{Proportion of target sequences with maximum two
#'   mismatches to at least one site outside the intended
#'   forward primer binding region in the input alignment.}
#'   \item{idOffTargetMatchFwd}{Names of all off-target matching sequences.}
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
#'   \item{offTargetMatchRev}{Proportion of target sequences with maximum two
#'   mismatches to at least one site outside the intended
#'   reverse primer binding region in the input alignment.}
#'   \item{idOffTargetMatchRev}{Names of all off-target matching sequences.}
#'
#' If the input assay contains probes, the following information
#' is also added:
#'
#'   \item{iupacSequencePr}{The probe sequence in IUPAC format.}
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
#'    \item{offTargetMatchPr}{Proportion of target sequences with maximum two
#'    mismatches to at least one site outside the intended
#'   probe binding region in the input alignment.}
#'   \item{idOffTargetMatchPr}{Names of all off-target matching sequences.}
#' }
#'
#' @export
#'
#' @references
#' Pages, H., Aboyoun, P., Gentleman R., and DebRoy S. (2020). Biostrings:
#' Efficient manipulation of biological strings. R package version
#' 2.57.2.
#'
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
#' #### RprimerAssay objects
#'
#' data("exampleRprimerAssay")
#' data("exampleRprimerAlignment")
#'
#' x <- exampleRprimerAssay[1:2, ]
#' target <- exampleRprimerAlignment
#'
#' checkMatch(x, target)
setGeneric("checkMatch", \(x, target) standardGeneric("checkMatch"))

# Methods ======================================================================

#' @describeIn checkMatch
#'
#' @export
setMethod("checkMatch", "RprimerOligo", \(x, target) {
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
setMethod("checkMatch", "RprimerAssay", \(x, target) {
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
#' target <- Biostrings::DNAStringSet(exampleRprimerAlignment)
#' .getSequenceNames(1:2, target)
.getSequenceNames <- function(x, target) {
    names(target) <- sub(" .*", "", names(target))
    lapply(x, \(i) names(target)[i])
}

#' @noRd
#'
#' @examples
#' data("exampleRprimerAlignment")
#' target <- exampleRprimerAlignment
#' x <- .maskRange(from = 5200, to = 5220, target, invert = TRUE)
.uniqueSequences <- function(x) {
    sequence <- unique(x)
    index <- match(unname(as.vector(x)), as.vector(sequence))
    name <- as.vector(names(x))
    names(name) <- index
    list(sequence = sequence, name = name)
}

#### continue from here ###############################################################

#x <- .maskRange(from = 5200, to = 5220, target, invert = TRUE)
#x <- .uniqueSequences(x)
#i <- .getMatchIndex(x$sequence, Biostrings::DNAStringSet(target), offTarget = FALSE)
# baka in unique sequences nanstans
# sen reduce...

########################################################################################

#' @noRd
#'
#' @examples
#' data("exampleRprimerAlignment")
#' data("exampleRprimerOligo")
#' x <- exampleRprimerOligo$sequence[[1]]
#' x <- Biostrings::DNAStringSet(x)
#' target <- Biostrings::DNAStringSet(exampleRprimerAlignment)
#' .getMatchIndex(x, target)
.getMatchIndex <- function(x, target, offTarget = FALSE) {
    if (!offTarget) {
        ## 0-3: number of mismatches we are looking for
        res <- lapply(seq(0, 3), \(i) {
            result <- Biostrings::vcountPDict(
                x, target,
                max.mismatch = i
            )
            which(colSums(result) == 0) ## Sequences who "fails"
        })
        res[-length(res)] <- lapply(seq(1, length(res) - 1), \(i) {
            setdiff(res[[i]], res[[i + 1]])
        })
        res[[length(res) + 1]] <- setdiff(seq_along(target), unlist(res))
        res[c(length(res), seq_len(length(res) - 1))]
    } else {
        result <- Biostrings::vcountPDict(
            x, target,
            max.mismatch = 2
        )
        list(which(colSums(result) > 0))
    }
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
.getMatchProportion <- function(x, target, offTarget = FALSE) {
    matching <- .getMatchIndex(x, target, offTarget)
    sequenceNames <- .getSequenceNames(matching, target)
    matching <- lapply(matching, \(x) length(x) / length(target))
    matching <- do.call("cbind.data.frame", matching)
    sequenceNames <- as.data.frame(t(do.call("cbind", list(sequenceNames))))
    matching <- cbind(matching, sequenceNames)
    if (!offTarget) {
        matching <- matching[c(1, 6, 2, 7, 3, 8, 4, 9, 5, 10)]
        names(matching) <- c(
            "perfectMatch", "idPerfectMatch",
            "oneMismatch", "idOneMismatch",
            "twoMismatches", "idTwoMismatches",
            "threeMismatches", "idThreeMismatches",
            "fourOrMoreMismatches", "idFourOrMoreMismatches"
        )
    } else {
        names(matching) <- c("offTargetMatch", "idOffTargetMatch")
    }
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
    onTarget <- lapply(seq_len(nrow(x)), \(i) {
        target <- .maskRange(
            x$start[[i]], x$end[[i]], target,
            invert = TRUE
        )
        check <- Biostrings::DNAStringSet(x$sequence[[i]])
        .getMatchProportion(check, target)
    })
    offTarget <- lapply(seq_len(nrow(x)), \(i) {
        target <- .maskRange(
            x$start[[i]], x$end[[i]], target,
            invert = FALSE
        )
        check <- Biostrings::DNAStringSet(x$sequence[[i]])
        .getMatchProportion(check, target, offTarget = TRUE)
    })
    onTarget <- do.call("rbind", onTarget)
    offTarget <- do.call("rbind", offTarget)
    cbind(
        "iupacSequence" = x$iupacSequence,
        onTarget, offTarget
    )
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
    rev$sequence <- lapply(rev$sequence, \(x) {
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
