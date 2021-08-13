# Helper functions to shiny app ================================================

filterOligos <- function(x,
                         fwdFrom,
                         fwdTo,
                         revFrom,
                         revTo,
                         prFrom,
                         prTo,
                         fwdIdentity,
                         fwdCoverage,
                         revIdentity,
                         revCoverage,
                         prIdentity,
                         prCoverage) {
    x <- as.data.frame(x)
    emptyRow <- makeEmptyRow(x)
    primers <- x[x$type == "primer", ]
    fwd <- primers$fwd &
        primers$start >= fwdFrom &
        primers$end <= fwdTo &
        primers$identity >= fwdIdentity &
        primers$coverage >= fwdCoverage
    rev <- primers$rev &
        primers$start >= revFrom &
        primers$end <= revTo &
        primers$identity >= revIdentity &
        primers$coverage >= revCoverage
    primers$fwd <- fwd
    primers$rev <- rev
    if (any(x$type == "probe")) {
        probes <- x[x$type == "probe", ]
        pr <- probes$start >= prFrom &
            probes$end <= prTo &
            probes$identity >= prIdentity &
            probes$coverage >= prCoverage
        probes <- probes[pr, ]
    } else {
        probes <- NULL
    }
    all <- rbind(primers, probes)
    all <- all[all$fwd | all$rev, ]
    all <- all[order(all$start), ]
    if (nrow(all) == 0L) {
        all <- emptyRow
    }
    RprimerOligo(all)
}

roundDbls <- function(x) {
    x <- as.data.frame(x)
    y <- lapply(seq_len(ncol(x)), function(i) {
        out <- if (is.double(x[, i])) {
            round(x[, i], 4)
        } else if (is.list(x[, i])) {
            if (is.double(x[, i][[1]])) {
                list(round(x[, i][[1]], 4))
            } else {
                list(x[, i][[1]])
            }
        } else {
            x[, i]
        }
    })
    names(y) <- names(x)
    data.frame(do.call("cbind", y))
}

removeListColumns <- function(x) {
    lists <- vapply(seq_len(ncol(x)), function(i) {
        is.list(x[, i])
    }, logical(1L))
    x[, !lists]
}


splitAssayToList <- function(x) {
    x <- as.data.frame(x)
    fwd <- x[, grepl("Fwd", names(x))]
    names(fwd) <- gsub("Fwd", "", names(fwd))
    rev <- x[, grepl("Rev", names(x))]
    names(rev) <- gsub("Rev", "", names(rev))
    if (any(grepl("Pr", names(x)))) {
        pr <- x[, grepl("Pr", names(x))]
        names(pr) <- gsub("Pr", "", names(pr))
        all <- list(fwd, rev, pr)
    } else {
        all <- list(fwd, rev)
    }
    all
}

makeListTable <- function(x) {
    x <- as.data.frame(x)
    list <- vapply(seq_len(ncol(x)), function(i) {
        is.list(x[, i])
    }, logical(1L))
    x <- x[, list]
    allCols <- lapply(seq_len(ncol(x)), function(i) {
        x[, i][[1]]
    })
    all <- data.frame(do.call("cbind", allCols))
    names(all) <- names(x)
    all$gcContent <- as.double(all$gcContent)
    all$tm <- as.double(all$tm)
    all$deltaG <- as.double(all$deltaG)
    all
}

makeEmptyRow <- function(x) {
    emptyRow <- x[1, , drop = FALSE]
    emptyRow[1, ] <- NA
    emptyRow$roiStart <- x$roiStart[[1]]
    emptyRow$roiEnd <- x$roiEnd[[1]]
    emptyRow$start <- 1
    emptyRow$end <- 1
    emptyRow$fwd <- emptyRow$rev <- TRUE
    emptyRow$type <- "primer"
    emptyRow$length <- NA
    emptyRow
}

assayOverviewTable <- function(x) {
    x <- as.data.frame(x)
    if (any(grepl("Pr", names(x)))) {
        x$iupacSequencePr <- ifelse(x$plusPr, x$iupacSequencePr, NA)
        x$iupacSequenceRcPr <- ifelse(x$minusPr, x$iupacSequenceRcPr, NA)
        x <- x[c(
            "start", "end", "length", "iupacSequenceFwd", "iupacSequenceRev",
            "iupacSequencePr", "iupacSequenceRcPr"
        )]
        names(x) <- c(
            "Start", "End", "Length", "Forward", "Reverse",
            "Probe, plus", "Probe, minus"
        )
    } else {
        x <- x[c(
            "start", "end", "length", "iupacSequenceFwd", "iupacSequenceRev"
        )]
        names(x) <- c(
            "Start", "End", "Length", "Forward", "Reverse"
        )
    }
    x
}
