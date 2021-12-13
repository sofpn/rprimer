## Helper functions to shiny app ===============================================

## UI functions ================================================================

switchPage <- function(i) {
    shiny::updateTabsetPanel(inputId = "wizard", selected = paste0("page", i))
}

previousButton <- function(i) {
    id <- paste0("page", i + 1, i)
    shiny::actionButton(
        id, "Previous",
        icon = shiny::icon("chevron-left")
    )
}

nextButton <- function(i) {
    id <- paste0("page", i - 1, i)
    shiny::actionButton(
        id,
        label = shiny::div(
            "Next", shiny::icon("chevron-right"),
            class = "btn btn-primary"
        )
    )
}

spinnerPlot <- function(id, ...) {
    shinycssloaders::withSpinner(
        shiny::plotOutput(id, ...),
        color = "grey"
    )
}

numericInputFrom <- function(x, id) {
    shiny::numericInput(id, "From",
        min = x$roiStart[[1]],
        max = x$roiEnd[[1]],
        value = x$roiStart[[1]]
    )
}

numericInputTo <- function(x, id) {
    shiny::numericInput(id, "To",
        min = x$roiStart[[1]],
        max = x$roiEnd[[1]],
        value = x$roiEnd[[1]]
    )
}

conservationInput <- function(x,
                              id,
                              type = "primer",
                              direction = "both",
                              variable = "identity") {
    x <- x[x$type == type, ]
    if (direction == "fwd") {
        x <- x[x$fwd, ]
    } else if (direction == "rev") {
        x <- x[x$rev, ]
    } else {
        x <- x[x$rev & x$fwd, ]
    }
    if (variable == "identity") {
        var <- x$identity
    } else if (variable == "coverage") {
        var <- x$coverage
    }
    minValue <- round(min(var, na.rm = TRUE), 4)
    maxValue <- round(max(var, na.rm = TRUE), 4)
    shiny::sliderInput(id,
        round = -4, step = 0.0001,
        paste("Minimum", variable),
        min = minValue,
        max = maxValue,
        value = minValue
    )
}

## Other functions =============================================================

arrangeMask <- function(from, to, aln) {
    alnLength <- ncol(aln)
    if (from == 0 && to == 0) {
        NULL
    } else {
        if (from < 1) from <- 1
        if (to < 1) to <- 1
        if (to > alnLength) to <- alnLength
        sort(from:to)
    }
}

filterOligos <- function(x,
                         fwdFrom = 0,
                         fwdTo = Inf,
                         revFrom = 0,
                         revTo = Inf,
                         prFrom = 0,
                         prTo = Inf,
                         fwdIdentity = 0,
                         fwdCoverage = 0,
                         revIdentity = 0,
                         revCoverage = 0,
                         prIdentity = 0,
                         prCoverage = 0) {
    x <- as.data.frame(x)
    emptyRow <- makeEmptyRow(x)
    if (any(x$type == "primer")) {
        primers <- x[x$type == "primer", ]
        if (any(primers$fwd)) {
            fwd <- primers$fwd &
                primers$start >= fwdFrom &
                primers$end <= fwdTo &
                primers$identity >= fwdIdentity &
                primers$coverage >= fwdCoverage
            primers$fwd <- fwd
        }
        if (any(primers$rev)) {
            rev <- primers$rev &
                primers$start >= revFrom &
                primers$end <= revTo &
                primers$identity >= revIdentity &
                primers$coverage >= revCoverage
            primers$rev <- rev
        }
    } else {
        primers <- NULL
    }
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
    y <- lapply(seq_len(ncol(x)), \(i) {
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
    lists <- vapply(seq_len(ncol(x)), \(i) {
        is.list(x[, i])
    }, logical(1L))
    x[, !lists]
}

splitAssayToList <- function(x) {
    x <- as.data.frame(x)
    fwd <- x[, grepl("Fwd", names(x))]
    fwd <- cbind(
        fwd,
        "roiStart" = x$roiStart, "roiEnd" = x$roiEnd
    )
    names(fwd) <- gsub("Fwd", "", names(fwd))
    rev <- x[, grepl("Rev", names(x))]
    rev <- cbind(
        rev,
        "roiStart" = x$roiStart, "roiEnd" = x$roiEnd
    )
    names(rev) <- gsub("Rev", "", names(rev))
    if (any(grepl("Pr", names(x)))) {
        pr <- x[, grepl("Pr", names(x))]
        pr <- cbind(
            pr,
            "roiStart" = x$roiStart, "roiEnd" = x$roiEnd
        )
        names(pr) <- gsub("Pr", "", names(pr))
        all <- list(fwd, rev, pr)
    } else {
        all <- list(fwd, rev)
    }
    all
}

convertToOligo <- function(x, rev = FALSE, type = "primer") {
    add <- data.frame(
        "type" = type, "fwd" = !rev, "rev" = rev, "score" = NA,
        "roiStart" = NA, "roiEnd" = NA, "iupacSequenceRc" = NA,
        "sequenceRc" = NA
    )
    x <- cbind(x, add)
    x$iupacSequenceRc <- makeRc(x$iupacSequence)
    x$sequenceRc <- makeRcList(x$sequence)
    if (rev) {
        old <- c("iupacSequence", "sequence", "iupacSequenceRc", "sequenceRc")
        new <- c("iupacSequenceRc", "sequenceRc", "iupacSequence", "sequence")
        names(x)[names(x) %in% old] <- new
    }
    x <- .beautifyOligos(x)
    RprimerOligo(x)
}

makeRc <- function(x) {
    x <- Biostrings::DNAString(x)
    x <- Biostrings::reverseComplement(x)
    as.character(x)
}

makeRcList <- function(x) {
    x <- unlist(x)
    x <- lapply(x, makeRc)
    x <- unlist(x)
    list(x)
}

makeListTable <- function(x) {
    x <- as.data.frame(x)
    list <- vapply(seq_len(ncol(x)), \(i) {
        is.list(x[, i])
    }, logical(1L))
    x <- x[, list]
    allCols <- lapply(seq_len(ncol(x)), \(i) {
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
