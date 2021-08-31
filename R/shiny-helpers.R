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
    shiny::numericInput(id, shiny::h5("From"),
                        min = x$roiStart[[1]],
                        max = x$roiEnd[[1]],
                        value = x$roiStart[[1]]
    )
}

numericInputTo <- function(x, id) {
    shiny::numericInput(id, shiny::h5("To"),
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
                       shiny::h5(paste("Minimum", variable)),
                       min = minValue,
                       max = maxValue,
                       value = minValue
    )
}


## Server functions ============================================================

displayDownloadHandlerTxt <- function(x, session) {
    shiny::renderUI({
        shiny::req(x)
        ns <- session$ns
        list(
            shiny::downloadLink(
                ns("downloadTxt"), "Download table as .txt"
            ),
            shiny::br()
        )
    })
}

displayDownloadHandlerFasta <- function(x, session) {
    shiny::renderUI({
        shiny::req(x)
        ns <- session$ns
        list(
            shiny::downloadLink(
                ns("downloadFasta"), "Download sequence(s) in fasta-format"
            ),
            shiny::br()
        )
    })
}

downloadHandlerTxt <- function(data, prefixName = NULL) {
    shiny::downloadHandler(
        filename <- function() {
            paste0(prefixName, "-", Sys.Date(), ".txt")
        },
        content <- function(file) {
            utils::write.table(
                as.data.frame(data), file,
                quote = FALSE, sep = "\t",
                row.names = FALSE
            )
        }
    )
}

downloadHandlerFasta <- function(data, prefixName = NULL) {
    shiny::downloadHandler(
        filename <- function() {
            paste0(prefixName, "-fasta-", Sys.Date(), ".txt")
        },
        content <- function(file) {
            data <- as(data, "DNAStringSet")
            Biostrings::writeXStringSet(data, file)
        }
    )
}

consensusDataTable <- function(data) {
    DT::renderDataTable(
        {
            shiny::req(is(data, "RprimerProfile"))
            x <- roundDbls(as.data.frame(data))
            names(x) <- c(
                "Position", "A", "C", "G", "T", "Other", "Gaps", "Majority",
                "Identity", "IUPAC", "Entropy", "Coverage"
            )
            x
        },
        options = list(
            info = FALSE,
            searching = FALSE, paging = TRUE,
            scrollX = TRUE, autoWidth = FALSE,
            ordering = FALSE, scrollY = "300"
        ),
        rownames = FALSE,
        selection = "none"
    )
}

oligoDataTable <- function(x, selection = "none", ordering = FALSE, ...) {
    DT::renderDataTable(
        {
            shiny::req(is(x, "RprimerOligo"))
            x <- roundDbls(removeListColumns(as.data.frame(x)))
            names(x) <- c(
                "Type", "Forward", "Reverse", "Start", "End", "Length",
                "IUPAC sequence", "IUPAC sequence, RC", "Identity",
                "Coverage", "Degeneracy", "GC content, mean",
                "GC content, range", "Tm, mean", "Tm, range", "Delta G, mean",
                "Delta G, range", "Design method", "Score", "ROI, start",
                "ROI, end"
            )
            if (is.na(x$Score[[1]])) {
                x <- x[!names(x) %in% "Score"]
            }
            x
        },
        options = list(
            info = FALSE,
            searching = FALSE, paging = FALSE,
            scrollX = TRUE, autoWidth = TRUE,
            ordering = ordering, ...
        ),
        rownames = FALSE,
        selection = selection
    )
}

assayDataTable <- function(x, selection = "none", ordering = FALSE, ...) {
    DT::renderDataTable(
        {
            shiny::req(!is.na(x$length[[1]]))
            x <- roundDbls(removeListColumns(as.data.frame(x)))
            if (any(grepl("Pr", names(x)))) {
                x$iupacSequencePr <- ifelse(x$plusPr, x$iupacSequencePr, NA)
                x$iupacSequenceRcPr <- ifelse(x$minusPr, x$iupacSequenceRcPr, NA)
            }
            names(x) <- if (any(grepl("Pr", names(x)))) {
                c(
                    "Start", "End", "Length", "Total degeneracy", "Score",
                    "Start, forward", "End, forward", "Length, forward",
                    "IUPAC sequence, forward", "Identity, forward",
                    "Coverage, forward", "Degeneracy, forward",
                    "GC content, mean, forward", "GC content, range, forward",
                    "Tm, mean, forward", "Tm, range, forward",
                    "Delta G, mean, forward", "Delta G, range, forward",
                    "Design method, forward",
                    "Start, reverse", "End, reverse", "Length, reverse",
                    "IUPAC sequence, reverse", "Identity, reverse",
                    "Coverage, reverse", "Degeneracy, reverse",
                    "GC content, mean, reverse", "GC content, range, reverse",
                    "Tm, mean, reverse", "Tm, range, reverse",
                    "Delta G, mean, reverse", "Delta G, range, reverse",
                    "Design method, reverse",
                    "Plus sense, probe", "Minus sense, probe",
                    "Start, probe", "End, probe", "Length, probe",
                    "IUPAC sequence, probe", "IUPAC sequence RC, probe",
                    "Identity, probe",
                    "Coverage, probe", "Degeneracy, probe",
                    "GC content, mean, probe", "GC content, range, probe",
                    "Tm, mean, probe", "Tm, range, probe",
                    "Delta G, mean, probe", "Delta G, range, probe",
                    "Design method, probe", "ROI, start", "ROI, end"
                )
            } else {
                c(
                    "Start", "End", "Length", "Total degeneracy", "Score",
                    "Start, forward", "End, forward", "Length, forward",
                    "IUPAC sequence, forward", "Identity, forward",
                    "Coverage, forward", "Degeneracy, forward",
                    "GC content, mean, forward", "GC content, range, forward",
                    "Tm, mean, forward", "Tm, range, forward",
                    "Delta G, mean, forward", "Delta G, range, forward",
                    "Design method, forward",
                    "Start, reverse", "End, reverse", "Length, reverse",
                    "IUPAC sequence, reverse", "Identity, reverse",
                    "Coverage, reverse", "Degeneracy, reverse",
                    "GC content, mean, reverse", "GC content, range, reverse",
                    "Tm, mean, reverse", "Tm, range, reverse",
                    "Delta G, mean, reverse", "Delta G, range, reverse",
                    "Design method, reverse",
                    "ROI, start", "ROI, end"
                )
            }
            x
        },
        options = list(
            info = FALSE,
            searching = FALSE, paging = TRUE,
            scrollX = TRUE, autoWidth = TRUE,
            ordering = ordering, ...
        ),
        rownames = FALSE,
        selection = selection
    )
}


assayOverviewTable <- function(x) {
    DT::renderDataTable(
        {
            shiny::req(x)
            if (is.na(x$length[[1]])) {
                NULL
            } else {
                x <- as.data.frame(x)
                if (any(grepl("Pr", names(x)))) {
                    x$iupacSequencePr <- ifelse(x$plusPr, x$iupacSequencePr, NA)
                    x$iupacSequenceRcPr <- ifelse(x$minusPr, x$iupacSequenceRcPr, NA)
                    x <- x[c(
                        "start", "end", "length", "iupacSequenceFwd",
                        "iupacSequenceRev",
                        "iupacSequencePr", "iupacSequenceRcPr"
                    )]
                    names(x) <- c(
                        "Start", "End", "Length", "Forward", "Reverse",
                        "Probe, plus", "Probe, minus"
                    )
                } else {
                    x <- x[c(
                        "start", "end", "length",
                        "iupacSequenceFwd", "iupacSequenceRev"
                    )]
                    names(x) <- c(
                        "Start", "End", "Length", "Forward", "Reverse"
                    )
                }
                x
            }
        },
        options = list(
            info = FALSE,
            searching = FALSE, paging = FALSE,
            scrollX = TRUE, autoWidth = FALSE,
            ordering = FALSE
        ),
        rownames = FALSE,
        selection = "none"
    )
}


oligoAllVariantTable <- function(x) {
    DT::renderDataTable(
        {
            shiny::req(x)
            x <- roundDbls(makeListTable(as.data.frame(x)))
            if (ncol(x) == 4) {
                names(x) <- c(
                    "Sequence", "GC content", "Tm", "Delta G"
                )
            } else {
                names(x) <- c(
                    "Sequence", "Sequence, RC", "GC content", "Tm", "Delta G"
                )
            }

            x
        },
        options = list(
            info = FALSE,
            searching = FALSE, paging = FALSE,
            scrollX = TRUE, autoWidth = FALSE,
            ordering = FALSE
        ),
        rownames = FALSE,
        selection = "none"
    )
}

oligoMatchTable <- function(x) {
    DT::renderDataTable(
        {
            shiny::req(x)
            x <- roundDbls(removeListColumns(as.data.frame(x)))
            names(x) <- c(
                "IUPAC sequence", "Perfect match", "1 mismatch", "2 mismatches",
                "3 mismatches", "4 or more mismatches"
            )
            x
        },
        options = list(
            info = FALSE,
            searching = FALSE, paging = FALSE,
            scrollX = TRUE, autoWidth = FALSE,
            ordering = FALSE
        ),
        rownames = FALSE,
        selection = "none"
    )
}

nucleotidePlot <- function(oligo, consensus, rc = FALSE) {
    shiny::renderPlot({
        shiny::req(oligo)
        shiny::req(consensus)
        from <- oligo$start
        to <- oligo$end
        bindingRegion <- consensus[
            consensus$position >= from & consensus$position <= to,
        ]
        plotData(bindingRegion, type = "nucleotide", rc)
    })
}

printMatchId <- function(x) {
    shiny::renderText({
        shiny::req(x)
        c(
            "<b>Perfect match</b><br>",
            x$idPerfectMatch[[1]],
            "<br><br><b>One mismatch</b><br>",
            x$idOneMismatch[[1]],
            "<br><br><b>Two mismatches</b><br>",
            x$idTwoMismatches[[1]],
            "<br><br><b>Three mismatches</b><br>",
            x$idThreeMismatches[[1]],
            "<br><br><b>Four or more mismatches</b><br>",
            x$idFourOrMoreMismatches[[1]],
            "<br><br><b>Off target match (< 3 mismatches)</b><br>",
            x$idOffTargetMatch[[1]]
        )
    })
}

## Other =======================================================================

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
    fwd <- cbind(
        fwd, "roiStart" = x$roiStart, "roiEnd" = x$roiEnd
    )
    names(fwd) <- gsub("Fwd", "", names(fwd))
    rev <- x[, grepl("Rev", names(x))]
    rev <- cbind(
        rev, "roiStart" = x$roiStart, "roiEnd" = x$roiEnd
    )
    names(rev) <- gsub("Rev", "", names(rev))
    if (any(grepl("Pr", names(x)))) {
        pr <- x[, grepl("Pr", names(x))]
        pr <- cbind(
            pr, "roiStart" = x$roiStart, "roiEnd" = x$roiEnd
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
