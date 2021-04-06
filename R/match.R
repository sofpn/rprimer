
.checkMatch <- function(x, target) {



    catchThreeEnd <- function(x) {

        vapply(x, function(y) {

            y <- unlist(strsplit(y, split = ""))

            y <- y[(length(y) - 7):length(y)] ## regulate the "size" of the 3' end here

            paste(y, collapse = "")

        }, character(1L), USE.NAMES = FALSE)

    }



    checkMatch <- function(x, target, max.mismatch = 4, threeEndsOnly = FALSE) {

        target <- Biostrings::DNAStringSet(target)

        target <- ShortRead::clean(target) ## Remove target sequences with wobble bases

        fwd <- unlist(x$sequenceFwd)

        rev <- unlist(x$sequenceRev)

        if (threeEndsOnly) {

            fwd <- catchThreeEnd(fwd)

            rev <- catchThreeEnd(rev)

        }

        fwd <- Biostrings::DNAStringSet(fwd)

        rev <- Biostrings::DNAStringSet(rev)

        rev <- Biostrings::reverseComplement(rev)

        matchFwd <- Biostrings::vcountPDict(fwd, target, max.mismatch = max.mismatch)

        matchRev <- Biostrings::vcountPDict(rev, target, max.mismatch = max.mismatch)

        isMatchFwd <- colSums(matchFwd) > 0

        isMatchRev <- colSums(matchRev) > 0

        isMatchBoth <- isMatchFwd & isMatchRev

        propMatchFwd <- sum(isMatchFwd) / length(isMatchFwd)

        propMatchRev <- sum(isMatchRev) / length(isMatchRev)

        propMatchBoth <- sum(isMatchBoth) / length(isMatchBoth)

        c("fwd" = propMatchFwd, "rev" = propMatchRev, "both" = propMatchBoth)

    }



    checkAssayMatch <- function(x,

                                target,

                                max.mismatch = 4,

                                threeEndsOnly = FALSE) {

        match <- lapply(seq_len(nrow(x)), function(i) {

            checkMatch(x[i, ], target, max.mismatch, threeEndsOnly)

        })

        match <- round(do.call("rbind", match), 2)

        rownames(match) <- paste0("assay_", seq_len(nrow(match)))

        match

    }



    checkAssayInclusitivty <- function(

        x, target, n.mismatch = 0:3, threeEndsOnly = FALSE

    ) {

        result <- lapply(n.mismatch, function(i) {

            checkAssayMatch(x, target, i, threeEndsOnly)

        })

        names(result) <- paste0("With_<=_", n.mismatch, "_mismatches")

        result

    }


}
