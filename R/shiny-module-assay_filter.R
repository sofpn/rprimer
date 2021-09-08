assayFilterUI <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h4("Filter assays (step 5/5)"),
        shiny::hr(),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::uiOutput(ns("assayFilter"))
            ),
            shiny::mainPanel(
                shiny::tabsetPanel(
                    id = ns("wizard"),
                    shiny::tabPanel(
                        title = "Plot",
                        shiny::br(),
                        spinnerPlot(
                            ns("assayPlot"),
                            width = "100%", height = 600
                        )
                    ),
                    shiny::tabPanel(
                        title = "Table",
                        shiny::br(),
                        shiny::uiOutput(ns("getDownloadLinkTxt")),
                        shiny::uiOutput(ns("getDownloadLinkFasta")),
                        shiny::br(),
                        shiny::h5("Select a row for more details"),
                        shiny::br(),
                        DT::dataTableOutput(ns("assayTable"))
                    ),
                    shiny::tabPanel(
                        title = "Selection",
                        shiny::br(),
                        shiny::h5("Assay information"),
                        shiny::br(),
                        shiny::uiOutput(ns("getDownloadLinkTxtSel")),
                        shiny::uiOutput(ns("getDownloadLinkFastaSel")),
                        shiny::br(),
                        DT::dataTableOutput(ns("overviewTable")),
                        shiny::br(),
                        shiny::br(),
                        spinnerPlot(ns("assayPositionPlot")),
                        shiny::br(),
                        shiny::h5("Amplicon sequence"),
                        shiny::hr(),
                        shiny::verbatimTextOutput(ns("ampliconSequence")),
                        shiny::br(),
                        shiny::h5("Oligo details"),
                        shiny::hr(),
                        shiny::uiOutput(ns("detailsTab"))
                    )
                )
            )
        ),
        shiny::hr()
    )
}


assayFilterServer <- function(id, alignment, consensus, oligo, allAssays) {
    shiny::moduleServer(id, function(input, output, session) {
        assay <- shiny::reactive({
            shiny::req(allAssays())
            x <- allAssays()
            x <- as.data.frame(x)
            emptyRow <- makeEmptyRow(x)
            x <- x[
                x$start >= input$assayRegionFrom & x$end <= input$assayRegionTo,
            ]
            x <- x[x$score <= input$maxAssayScore, ]
            if (nrow(x) == 0L) {
                x <- emptyRow
            }
            RprimerAssay(x)
        })

        output$assayFilter <- shiny::renderUI({
            ns <- session$ns

            list(
                numericInputFrom(allAssays(), ns("assayRegionFrom")),
                numericInputTo(allAssays(), ns("assayRegionTo")),
                sliderInput(ns("maxAssayScore"),
                    round = -4, step = 0.0001,
                    shiny::h5("Maximum score (lower is better)"),
                    min = round(min(allAssays()$score, na.rm = TRUE), 4),
                    max = round(max(allAssays()$score, na.rm = TRUE), 4),
                    value = round(max(allAssays()$score, na.rm = TRUE), 4)
                )
            )
        })

        output$assayPlot <- shiny::renderPlot({
            shiny::req(is(assay(), "RprimerAssay"))
            plotData(assay())
        })


        output$getDownloadLinkTxt <- shiny::renderUI({
            shiny::req(is(assay(), "RprimerAssay"))
            ns <- session$ns
            list(
                shiny::downloadLink(
                    ns("downloadTxt"), "Download table as .txt"
                ),
                shiny::br()
            )
        })

        output$getDownloadLinkFasta <- shiny::renderUI({
            shiny::req(is(assay(), "RprimerAssay"))
            ns <- session$ns
            list(
                shiny::downloadLink(
                    ns("downloadFasta"), "Download sequence(s) in fasta-format"
                ),
                shiny::br()
            )
        })

        output$downloadTxt <- shiny::downloadHandler(
            filename <- function() {
                paste0("assay", "-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                utils::write.table(
                    as.data.frame(assay()), file,
                    quote = FALSE, sep = "\t",
                    row.names = FALSE
                )
            }
        )

        output$downloadFasta <- shiny::downloadHandler(
            filename <- function() {
                paste0("assay", "-fasta-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                data <- as(assay(), "DNAStringSet")
                Biostrings::writeXStringSet(data, file)
            }
        )

        output$assayTable <- DT::renderDataTable(
            {
                shiny::req(!is.na(assay()$length[[1]]))
                x <- roundDbls(removeListColumns(as.data.frame(assay())))
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
                ordering = TRUE, scrollY = "300"
            ),
            rownames = FALSE,
            selection = list(mode = "single")
        )

        shiny::observeEvent(input$assayTable_rows_selected, {
            Sys.sleep(1)
            shiny::updateTabsetPanel(session, "wizard", selected = "Selection")
        })

        selectedAssay <- shiny::reactive({
            shiny::req(assay())
            if (!is.null(input$assayTable_rows_selected)) {
                assay()[input$assayTable_rows_selected, ]
            } else {
                NULL
            }
        })

        assayList <- shiny::reactive({
            shiny::req(selectedAssay())
            splitAssayToList(selectedAssay())
        })

        assayMatch <- shiny::reactive({
            shiny::req(selectedAssay())
            shiny::req(alignment())
            if (is.na(selectedAssay()$length[[1]])) {
                NULL
            } else {
                checkMatch(selectedAssay(), alignment())
            }
        })

        assayMatchList <- shiny::reactive({
            shiny::req(assayMatch())
            splitAssayToList(assayMatch())
        })

        output$ampliconSequence <- shiny::renderPrint({
            shiny::req(selectedAssay())
            from <- selectedAssay()$start
            to <- selectedAssay()$end
            sequence <- consensus()$iupac[
                consensus()$position >= from & consensus()$position <= to]
            paste(sequence, collapse = "")
        })

        output$getDownloadLinkTxtSel <- shiny::renderUI({
            shiny::req(is(selectedAssay(), "RprimerAssay"))
            ns <- session$ns
            list(
                shiny::downloadLink(
                    ns("downloadTxtSel"), "Download table as .txt"
                ),
                shiny::br()
            )
        })

        output$getDownloadLinkFastaSel <- shiny::renderUI({
            shiny::req(is(selectedAssay(), "RprimerAssay"))
            ns <- session$ns
            list(
                shiny::downloadLink(
                    ns("downloadFastaSel"), "Download sequence(s) in fasta-format"
                ),
                shiny::br()
            )
        })

        output$downloadTxtSel <- shiny::downloadHandler(
            filename <- function() {
                paste0("assay-selection", "-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                utils::write.table(
                    as.data.frame(selectedAssay()), file,
                    quote = FALSE, sep = "\t",
                    row.names = FALSE
                )
            }
        )

        output$downloadFastaSel <- shiny::downloadHandler(
            filename <- function() {
                paste0("assay-selection", "-fasta-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                data <- as(selectedAssay(), "DNAStringSet")
                Biostrings::writeXStringSet(data, file)
            }
        )

        output$overviewTable <- DT::renderDataTable(
            {
                shiny::req(is(selectedAssay(), "RprimerAssay"))
                if (is.na(selectedAssay()$length[[1]])) {
                    NULL
                } else {
                    x <- as.data.frame(selectedAssay())
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

        output$assayPositionPlot <- shiny::renderPlot({
            req(selectedAssay())
            req(consensus())
            plotData(
                consensus(),
                highlight = c(selectedAssay()$start, selectedAssay()$end)
            )
        })

        output$detailsTab <- shiny::renderUI({
            ns <- session$ns
            tabs <- list(
                shiny::tabPanel(
                    title = "Forward",
                    br(),
                    shiny::tabsetPanel(type = "pills",
                        shiny::tabPanel(
                            title = "Oligo information",
                            shiny::br(),
                            shiny::h5("Overview"),
                            shiny::hr(),
                            DT::dataTableOutput(ns("overviewTableFwd")),
                            shiny::br(),
                            shiny::h5("All sequence variants"),
                            shiny::hr(),
                            DT::dataTableOutput(ns("allVariantTableFwd")),
                            shiny::br(),
                            shiny::h5("Nucleotide distribution in target alignment"),
                            shiny::hr(),
                            shiny::br(),
                            shiny::column(
                                width = 12, align = "center",
                                spinnerPlot(
                                    ns("bindingRegionPlotFwd"),
                                    width = "75%"
                                )
                            )
                        ),
                        shiny::tabPanel(
                            title = "Match details",
                            shiny::br(),
                            shiny::h5("Proportion of matching sequences"),
                            shiny::br(),
                            DT::dataTableOutput(ns("matchTableFwd")),
                            shiny::br(),
                            shiny::column(
                                width = 12, align = "center",
                                spinnerPlot(ns("matchPlotFwd"), width = "75%")
                            ),
                            shiny::br(),
                            shiny::h5("Sequence names"),
                            shiny::hr(),
                            shiny::htmlOutput(ns("matchIdFwd"))
                        )
                    )
                ),
                shiny::tabPanel(
                    title = "Reverse",
                    br(),
                    shiny::tabsetPanel(type = "pills",
                        shiny::tabPanel(
                            title = "Oligo information",
                            shiny::br(),
                            shiny::h5("Overview"),
                            shiny::hr(),
                            DT::dataTableOutput(ns("overviewTableRev")),
                            shiny::br(),
                            shiny::h5("All sequence variants"),
                            shiny::hr(),
                            DT::dataTableOutput(ns("allVariantTableRev")),
                            shiny::br(),
                            shiny::h5("Nucleotide distribution in target alignment"),
                            shiny::hr(),
                            shiny::br(),
                            shiny::column(
                                width = 12, align = "center",
                                spinnerPlot(
                                    ns("bindingRegionPlotRev"),
                                    width = "75%"
                                )
                            )
                        ),
                        shiny::tabPanel(
                            title = "Match details",
                            shiny::br(),
                            shiny::h5("Proportion of matching sequences"),
                            shiny::br(),
                            DT::dataTableOutput(ns("matchTableRev")),
                            shiny::br(),
                            shiny::column(
                                width = 12, align = "center",
                                spinnerPlot(ns("matchPlotRev"), width = "75%")
                            ),
                            shiny::br(),
                            shiny::h5("Sequence names"),
                            shiny::hr(),
                            shiny::htmlOutput(ns("matchIdRev"))
                        )
                    )
                ),
                shiny::tabPanel(
                    title = "Probe",
                    br(),
                    shiny::tabsetPanel(type = "pills",
                        shiny::tabPanel(
                            title = "Oligo information",
                            shiny::br(),
                            shiny::h5("Overview"),
                            shiny::hr(),
                            DT::dataTableOutput(ns("overviewTablePr")),
                            shiny::br(),
                            shiny::h5("All sequence variants"),
                            shiny::hr(),
                            DT::dataTableOutput(ns("allVariantTablePr")),
                            shiny::br(),
                            shiny::h5("Nucleotide distribution in target alignment"),
                            shiny::hr(),
                            shiny::br(),
                            shiny::column(
                                width = 12, align = "center",
                                spinnerPlot(ns("bindingRegionPlotPr"), width = "75%")
                            )
                        ),
                        shiny::tabPanel(
                            title = "Match details",
                            shiny::br(),
                            shiny::h5("Proportion of matching sequences"),
                            shiny::br(),
                            DT::dataTableOutput(ns("matchTablePr")),
                            shiny::br(),
                            shiny::column(
                                width = 12, align = "center",
                                spinnerPlot(ns("matchPlotPr"), width = "75%")
                            ),
                            shiny::br(),
                            shiny::h5("Sequence names"),
                            shiny::hr(),
                            shiny::htmlOutput(ns("matchIdPr"))
                        )
                    )
                )
            )

            if (length(assayList()) == 3) {
                do.call(shiny::tabsetPanel, tabs)
            } else {
                do.call(shiny::tabsetPanel, tabs[seq_len(2)])
            }
        })

        oligoFwd <- reactive(convertToOligo(assayList()[[1]]))

        oligoRev <- reactive(convertToOligo(assayList()[[2]], rev = TRUE))

        oligoPr <- reactive(
            if (length(assayList()) == 3L) {
                x <- assayList()[[3]]
                convertToOligo(x, type = "probe")
            } else {
                NULL
            }
        )

        matchFwd <- shiny::reactive({
            shiny::req(oligoFwd())
            shiny::req(alignment())
            checkMatch(oligoFwd(), alignment())
        })

        output$overviewTableFwd <- DT::renderDataTable(
            {
                shiny::req(is(oligoFwd(), "RprimerOligo"))
                x <- roundDbls(removeListColumns(as.data.frame(oligoFwd())))
                names(x) <- c(
                    "Type", "Forward", "Reverse", "Start", "End", "Length",
                    "IUPAC sequence", "IUPAC sequence, RC", "Identity",
                    "Coverage", "Degeneracy", "GC content, mean",
                    "GC content, range", "Tm, mean", "Tm, range", "Delta G, mean",
                    "Delta G, range", "Design method", "Score", "ROI, start",
                    "ROI, end"
                )
                x <- x[!names(x) %in% c(
                    "Type", "Forward", "Reverse", "Score", "IUPAC sequence, RC",
                    "ROI, start", "ROI, end"
                )]
                x
            },
            options = list(
                info = FALSE,
                searching = FALSE, paging = FALSE,
                scrollX = TRUE, autoWidth = TRUE,
                ordering = FALSE
            ),
            rownames = FALSE,
            selection = "none"
        )

        output$allVariantTableFwd <- DT::renderDataTable(
            {
                shiny::req(is(oligoFwd(), "RprimerOligo"))
                x <- roundDbls(makeListTable(as.data.frame(oligoFwd())))
                names(x) <- c(
                    "Sequence", "Sequence, RC", "GC content", "Tm", "Delta G"
                )
                x <- x[!names(x) %in% c("Sequence, RC")]
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

        output$bindingRegionPlotFwd <- shiny::renderPlot({
            shiny::req(is(oligoFwd(), "RprimerOligo"))
            shiny::req(consensus())
            from <- oligoFwd()$start
            to <- oligoFwd()$end
            bindingRegion <- consensus()[
                consensus()$position >= from & consensus()$position <= to,
            ]
            plotData(bindingRegion, type = "nucleotide")
        })

        output$matchTableFwd <- DT::renderDataTable(
            {
                shiny::req(matchFwd())
                x <- roundDbls(removeListColumns(as.data.frame(matchFwd())))
                names(x) <- c(
                    "IUPAC sequence", "Perfect match", "1 mismatch", "2 mismatches",
                    "3 mismatches", "4 or more mismatches",
                    "Off target match (<3 mismatches)"
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

        output$matchPlotFwd <- shiny::renderPlot({
            shiny::req(matchFwd())
            plotData(matchFwd())
        })

        output$matchIdFwd <- shiny::renderText({
            shiny::req(matchFwd())
            x <- matchFwd()
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

        matchRev <- shiny::reactive({
            shiny::req(oligoRev())
            shiny::req(alignment())
            checkMatch(oligoRev(), alignment())
        })

        output$overviewTableRev <- DT::renderDataTable(
            {
                shiny::req(is(oligoRev(), "RprimerOligo"))
                x <- roundDbls(removeListColumns(as.data.frame(oligoRev())))
                names(x) <- c(
                    "Type", "Forward", "Reverse", "Start", "End", "Length",
                    "IUPAC sequence, RC", "IUPAC sequence", "Identity",
                    "Coverage", "Degeneracy", "GC content, mean",
                    "GC content, range", "Tm, mean", "Tm, range", "Delta G, mean",
                    "Delta G, range", "Design method", "Score", "ROI, start",
                    "ROI, end"
                )
                x <- x[!names(x) %in% c(
                    "Type", "Forward", "Reverse", "Score", "IUPAC sequence, RC",
                    "ROI, start", "ROI, end"
                )]
                x
            },
            options = list(
                info = FALSE,
                searching = FALSE, paging = FALSE,
                scrollX = TRUE, autoWidth = TRUE,
                ordering = FALSE
            ),
            rownames = FALSE,
            selection = "none"
        )

        output$allVariantTableRev <- DT::renderDataTable(
            {
                shiny::req(is(oligoRev(), "RprimerOligo"))
                x <- roundDbls(makeListTable(as.data.frame(oligoRev())))
                names(x) <- c(
                    "Sequence, RC", "Sequence", "GC content", "Tm", "Delta G"
                )
                x <- x[!names(x) %in% c("Sequence, RC")]
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

        output$bindingRegionPlotRev <- shiny::renderPlot({
            shiny::req(is(oligoRev(), "RprimerOligo"))
            shiny::req(consensus())
            from <- oligoRev()$start
            to <- oligoRev()$end
            bindingRegion <- consensus()[
                consensus()$position >= from & consensus()$position <= to,
            ]
            plotData(bindingRegion, type = "nucleotide", rc = TRUE)
        })

        output$matchTableRev <- DT::renderDataTable(
            {
                shiny::req(matchRev())
                x <- roundDbls(removeListColumns(as.data.frame(matchRev())))
                names(x) <- c(
                    "IUPAC sequence", "Perfect match", "1 mismatch", "2 mismatches",
                    "3 mismatches", "4 or more mismatches",
                    "Off target match (<3 mismatches)"
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

        output$matchPlotRev <- shiny::renderPlot({
            shiny::req(matchRev())
            plotData(matchRev())
        })

        output$matchIdRev <- shiny::renderText({
            shiny::req(matchRev())
            x <- matchRev()
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

        matchPr <- shiny::reactive({
            shiny::req(oligoPr())
            shiny::req(alignment())
            checkMatch(oligoPr(), alignment())
        })

        output$overviewTablePr <- DT::renderDataTable(
            {
                shiny::req(is(oligoPr(), "RprimerOligo"))
                x <- roundDbls(removeListColumns(as.data.frame(oligoPr())))
                names(x) <- c(
                    "Type", "Plus", "Minus", "Start", "End", "Length",
                    "IUPAC sequence, plus", "IUPAC sequence, minus", "Identity",
                    "Coverage", "Degeneracy", "GC content, mean",
                    "GC content, range", "Tm, mean", "Tm, range", "Delta G, mean",
                    "Delta G, range", "Design method", "Score", "ROI, start",
                    "ROI, end"
                )
                x <- x[!names(x) %in% c(
                    "Type", "Score", "ROI, start", "ROI, end"
                )]
                x
            },
            options = list(
                info = FALSE,
                searching = FALSE, paging = FALSE,
                scrollX = TRUE, autoWidth = TRUE,
                ordering = FALSE
            ),
            rownames = FALSE,
            selection = "none"
        )

        output$allVariantTablePr <- DT::renderDataTable(
            {
                shiny::req(is(oligoPr(), "RprimerOligo"))
                x <- roundDbls(makeListTable(as.data.frame(oligoPr())))
                names(x) <- c(
                    "Sequence, plus", "Sequence, minus", "GC content", "Tm", "Delta G"
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

        output$bindingRegionPlotPr <- shiny::renderPlot({
            shiny::req(is(oligoPr(), "RprimerOligo"))
            shiny::req(consensus())
            from <- oligoPr()$start
            to <- oligoPr()$end
            bindingRegion <- consensus()[
                consensus()$position >= from & consensus()$position <= to,
            ]
            plotData(bindingRegion, type = "nucleotide")
        })

        output$matchTablePr <- DT::renderDataTable(
            {
                shiny::req(matchPr())
                x <- roundDbls(removeListColumns(as.data.frame(matchPr())))
                names(x) <- c(
                    "IUPAC sequence", "Perfect match", "1 mismatch", "2 mismatches",
                    "3 mismatches", "4 or more mismatches",
                    "Off target match (<3 mismatches)"
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

        output$matchPlotPr <- shiny::renderPlot({
            shiny::req(matchPr())
            plotData(matchPr())
        })

        output$matchIdPr <- shiny::renderText({
            shiny::req(matchPr())
            x <- matchPr()
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

        list(data = shiny::reactive(assay()))
    })
}
