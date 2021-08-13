assayFilterUI <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::titlePanel("Step 5/5: Filter assays"),
        shiny::br(),
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
                        shinycssloaders::withSpinner(
                            shiny::plotOutput(
                                ns("assayPlot"),
                                width = "100%",
                                height = 600
                            ),
                            color = "grey"
                        )
                    ),
                    shiny::tabPanel(
                        title = "Table",
                        shiny::br(),
                        shiny::uiOutput(ns("getData")),
                        shiny::br(),
                        shiny::h5("Select a row for more details"),
                        shiny::br(),
                        DT::dataTableOutput(ns("assayTable"))
                    ),
                    shiny::tabPanel(
                        title = "Selection",
                        shiny::br(),
                        width = 12,
                        shiny::tabsetPanel(
                            shiny::tabPanel(
                                title = "Overview",
                                shiny::br(),
                                shiny::uiOutput(ns("getSelectionData")),
                                shiny::br(),
                                shiny::h5("General assay information"),
                                shiny::hr(),
                                DT::dataTableOutput(ns("assayTableSelection")),
                                shiny::br(),
                                shiny::h5("Details"),
                                shiny::hr(),
                                shiny::br(),
                                shiny::uiOutput(ns("detailsTab"))
                            ),
                            shiny::tabPanel(
                                title = "Match details",
                                shiny::br(),
                                shiny::h5(
                                    "Proportion of matching sequences within the
                                                  intended target binding region in the input alignment"
                                ),
                                shiny::hr(),
                                shiny::br(),
                                shiny::column(
                                    width = 12, align = "center",
                                    shinycssloaders::withSpinner(
                                        shiny::plotOutput(
                                            ns("matchPlot"),
                                            width = "100%"
                                        ),
                                        color = "grey"
                                    )
                                ),
                                shiny::h5("Details"),
                                shiny::hr(),
                                shiny::br(),
                                shiny::uiOutput(ns("matchDetailsTab"))
                            )
                        )
                    )
                )
            )
        )
    )
}

assayFilterServer <- function(id, alignment, consensus, assays) {
    shiny::moduleServer(id, function(input, output, session) {
        assay <- shiny::reactive({
            shiny::req(is(assays(), "RprimerAssay"))
            x <- assays()
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
                numericInput(ns("assayRegionFrom"), shiny::h5("From"),
                    min = assays()$roiStart[[1]],
                    max = assays()$roiEnd[[1]],
                    value = assays()$roiStart[[1]]
                ),
                numericInput(ns("assayRegionTo"), shiny::h5("To"),
                    min = assays()$roiStart[[1]],
                    max = assays()$roiEnd[[1]],
                    value = assays()$roiEnd[[1]]
                ),
                sliderInput(ns("maxAssayScore"),
                    shiny::h5("Maximum score (lower is better)"),
                    min = min(assays()$score, na.rm = TRUE),
                    max = max(assays()$score, na.rm = TRUE),
                    value = max(assays()$score, na.rm = TRUE)
                )
            )
        })

        output$assayPlot <- shiny::renderPlot({
            shiny::req(is(assay(), "RprimerAssay"))
            plotData(assay())
        })

        shiny::observeEvent(input$assayTable_rows_selected, {
            Sys.sleep(1)
            shiny::updateTabsetPanel(session, "wizard", selected = "Selection")
        })

        selectedAssay <- shiny::reactive({
            shiny::req(is(assay(), "RprimerAssay"))
            if (!is.null(input$assayTable_rows_selected)) {
                assay()[input$assayTable_rows_selected, ]
            } else {
                NULL
            }
        })

        selectedAssayList <- shiny::reactive({
            shiny::req(selectedAssay())
            splitAssayToList(selectedAssay())
        })

        selectedAssayMatch <- shiny::reactive({
            shiny::req(selectedAssay())
            shiny::req(alignment())
            if (is.na(selectedAssay()$length[[1]])) {
                NULL
            } else {
                checkMatch(selectedAssay(), alignment())
            }
        })

        selectedAssayMatchList <- shiny::reactive({
            shiny::req(selectedAssayMatch())
            splitAssayToList(selectedAssayMatch())
        })

        output$downloadTable <- shiny::downloadHandler(
            filename <- function() {
                paste0("assays-filtered-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                write.table(
                    as.data.frame(assaySelection()), file,
                    quote = FALSE, sep = "\t",
                    row.names = FALSE
                )
            }
        )

        output$getData <- shiny::renderUI({
            shiny::req(!is.na(assay()$length[[1]]))
            ns <- session$ns
            list(
                shiny::downloadLink(
                    ns("downloadTable"), "Download table as .txt"
                ),
                shiny::br(),
                shiny::downloadLink(
                    ns("downloadFasta"),
                    "Download assays in fasta-format"
                )
            )
        })

        output$downloadFasta <- shiny::downloadHandler(
            filename <- function() {
                paste0("assays-filtered-fasta", Sys.Date(), ".txt")
            },
            content <- function(file) {
                x <- as(assaySelection(), "DNAStringSet")
                Biostrings::writeXStringSet(x, file)
            }
        )


        output$assayTable <- DT::renderDataTable(
            {
                shiny::req(!is.na(assay()$length[[1]]))
                x <- assay()
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
                ordering = TRUE, scrollY = "300"
            ),
            rownames = FALSE,
            selection = list(mode = "single")
        )

        output$detailsTab <- shiny::renderUI({
            ns <- session$ns
            tabs <- list(
                shiny::tabPanel(
                    title = "Forward",
                    shiny::br(),
                    DT::dataTableOutput(ns("tableFwd")),
                    shiny::br(),
                    shiny::br(),
                    shiny::h5("All sequence variants"),
                    shiny::hr(),
                    shiny::br(),
                    DT::dataTableOutput(ns("tableFwdAll")),
                    shiny::br(),
                    shiny::br(),
                    shiny::h5("Nucleotide distribution in target alignment"),
                    shiny::hr(),
                    shiny::br(),
                    shiny::column(
                        width = 12, align = "center",
                        shinycssloaders::withSpinner(
                            shiny::plotOutput(
                                ns("ntPlotFwd"),
                                width = "75%"
                            ),
                            color = "grey"
                        )
                    ),
                    shiny::br(),
                    shiny::br()
                ),
                shiny::tabPanel(
                    title = "Reverse",
                    shiny::br(),
                    DT::dataTableOutput(ns("tableRev")),
                    shiny::br(),
                    shiny::br(),
                    shiny::h5("All sequence variants"),
                    shiny::hr(),
                    shiny::br(),
                    DT::dataTableOutput(ns("tableRevAll")),
                    shiny::br(),
                    shiny::br(),
                    shiny::h5("Nucleotide distribution in target alignment"),
                    shiny::hr(),
                    shiny::br(),
                    shiny::column(
                        width = 12, align = "center",
                        shinycssloaders::withSpinner(
                            shiny::plotOutput(
                                ns("ntPlotRev"),
                                width = "75%"
                            ),
                            color = "grey"
                        )
                    )
                ),
                shiny::tabPanel(
                    title = "Probe",
                    shiny::br(),
                    DT::dataTableOutput(ns("tablePr")),
                    shiny::br(),
                    shiny::br(),
                    shiny::h5("All sequence variants"),
                    shiny::hr(),
                    shiny::br(),
                    DT::dataTableOutput(ns("tablePrAll")),
                    shiny::br(),
                    shiny::br(),
                    shiny::h5("Nucleotide distribution in target alignment"),
                    shiny::hr(),
                    shiny::br(),
                    shiny::column(
                        width = 12, align = "center",
                        shinycssloaders::withSpinner(
                            shiny::plotOutput(
                                ns("ntPlotPr"),
                                width = "75%"
                            ),
                            color = "grey"
                        )
                    ),
                    shiny::br(),
                    shiny::br()
                )
            )

            if (length(selectedAssayList()) == 3) {
                do.call(shiny::tabsetPanel, tabs)
            } else {
                do.call(shiny::tabsetPanel, tabs[1:2])
            }
        })

        output$matchDetailsTab <- shiny::renderUI({
            ns <- session$ns

            tabs <- list(
                shiny::tabPanel(
                    title = "Forward",
                    shiny::br(),
                    DT::dataTableOutput(
                        ns("tableFwdMatch")
                    ),
                    shiny::br(),
                    shiny::h5("Target sequence names"),
                    shiny::hr(),
                    shiny::htmlOutput(ns("perfectMatchFwd")),
                    shiny::br(),
                    shiny::br(),
                    shiny::htmlOutput(ns("oneMismatchFwd")),
                    shiny::br(),
                    shiny::br(),
                    shiny::htmlOutput(ns("twoMismatchesFwd")),
                    shiny::br(),
                    shiny::br(),
                    shiny::htmlOutput(ns("threeMismatchesFwd")),
                    shiny::br(),
                    shiny::br(),
                    shiny::htmlOutput(ns("fourOrMoreMismatchesFwd")),
                    shiny::br(),
                    shiny::br()
                ),
                shiny::tabPanel(
                    title = "Reverse",
                    shiny::br(),
                    DT::dataTableOutput(
                        ns("tableRevMatch")
                    ),
                    shiny::br(),
                    shiny::h5("Target sequence names"),
                    shiny::hr(),
                    shiny::htmlOutput(ns("perfectMatchRev")),
                    shiny::br(),
                    shiny::br(),
                    shiny::htmlOutput(ns("oneMismatchRev")),
                    shiny::br(),
                    shiny::br(),
                    shiny::htmlOutput(ns("twoMismatchesRev")),
                    shiny::br(),
                    shiny::br(),
                    shiny::htmlOutput(ns("threeMismatchesRev")),
                    shiny::br(),
                    shiny::br(),
                    shiny::htmlOutput(ns("fourOrMoreMismatchesRev")),
                    shiny::br(),
                    shiny::br()
                ),
                shiny::tabPanel(
                    title = "Probe",
                    shiny::br(),
                    DT::dataTableOutput(
                        ns("tablePrMatch")
                    ),
                    shiny::br(),
                    shiny::h5("Target sequence names"),
                    shiny::hr(),
                    shiny::htmlOutput(ns("perfectMatchPr")),
                    shiny::br(),
                    shiny::br(),
                    shiny::htmlOutput(ns("oneMismatchPr")),
                    shiny::br(),
                    shiny::br(),
                    shiny::htmlOutput(ns("twoMismatchesPr")),
                    shiny::br(),
                    shiny::br(),
                    shiny::htmlOutput(ns("threeMismatchesPr")),
                    shiny::br(),
                    shiny::br(),
                    shiny::htmlOutput(ns("fourOrMoreMismatchesPr")),
                    shiny::br(),
                    shiny::br()
                )
            )

            if (length(selectedAssayList()) == 3) {
                do.call(shiny::tabsetPanel, tabs)
            } else {
                do.call(shiny::tabsetPanel, tabs[1:2])
            }
        })





        output$assayTableSelection <- DT::renderDataTable(
            {
                shiny::req(selectedAssay())
                if (is.na(selectedAssay()$length[[1]])) {
                    NULL
                } else {
                    assayOverviewTable(selectedAssay())
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

        output$ntPlotFwd <- shiny::renderPlot({
            shiny::req(selectedAssay())
            if (is.na(selectedAssay()$length[[1]])) {
                NULL
            } else {
                from <- selectedAssay()$startFwd
                to <- selectedAssay()$endFwd
                plotData(consensus()[
                    consensus()$position >= from & consensus()$position <= to,
                ], type = "nucleotide")
            }
        })

        output$ntPlotRev <- shiny::renderPlot({
            shiny::req(selectedAssay())
            if (is.na(selectedAssay()$length[[1]])) {
                NULL
            } else {
                from <- selectedAssay()$startRev
                to <- selectedAssay()$endRev
                plotData(consensus()[
                    consensus()$position >= from & consensus()$position <= to,
                ], type = "nucleotide", rc = TRUE)
            }
        })

        output$ntPlotPr <- shiny::renderPlot({
            shiny::req(selectedAssay())
            shiny::req(grepl("Pr", names(selectedAssay())))
            if (is.na(selectedAssay()$length[[1]])) {
                NULL
            } else {
                from <- selectedAssay()$startPr
                to <- selectedAssay()$endPr
                rc <- ifelse(selectedAssay()$plusPr, FALSE, TRUE)
                plotData(consensus()[
                    consensus()$position >= from & consensus()$position <= to,
                ], type = "nucleotide", rc = rc)
            }
        })

        output$getSelectionData <- shiny::renderUI({
            shiny::req(selectedAssay())
            ns <- session$ns
            list(
                shiny::downloadLink(
                    ns("downloadSelectionTable"),
                    "Download assay information as .txt"
                ),
                shiny::br(),
                shiny::downloadLink(
                    ns("downloadSelectionFasta"),
                    "Download assay in fasta-format"
                )
            )
        })

        output$downloadSelectionTable <- shiny::downloadHandler(
            filename <- function() {
                paste0("assay-filtered-selection-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                write.table(
                    as.data.frame(selectedAssay()), file,
                    quote = FALSE, sep = "\t",
                    row.names = FALSE
                )
            }
        )

        output$downloadSelectionFasta <- shiny::downloadHandler(
            filename <- function() {
                paste0("assay-filtered-selection-fasta-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                x <- selectedAssay()
                x <- as(selectedAssay(), "DNAStringSet")
                Biostrings::writeXStringSet(x, file)
            }
        )

        output$tableFwd <- DT::renderDataTable(
            {
                x <- roundDbls(removeListColumns(selectedAssayList()[[1]]))
                names(x) <- c(
                    "Start", "End", "Length", "IUPAC sequence", "Identity",
                    "Coverage", "Degeneracy", "GC, mean", "GC, range", "Tm, mean",
                    "Tm, range", "Delta G, mean", "Delta G, range", "Method"
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


        output$tableFwdAll <- DT::renderDataTable(
            {
                shiny::req(selectedAssayList())
                x <- roundDbls(makeListTable(as.data.frame(selectedAssayList()[[1]])))
                names(x) <- c(
                    "Sequence", "GC content", "Tm", "Delta G"
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

        output$tableRev <- DT::renderDataTable(
            {
                shiny::req(selectedAssayList())
                x <- roundDbls(removeListColumns(selectedAssayList()[[2]]))
                names(x) <- c(
                    "Start", "End", "Length", "IUPAC sequence", "Identity",
                    "Coverage", "Degeneracy", "GC, mean", "GC, range", "Tm, mean",
                    "Tm, range", "Delta G, mean", "Delta G, range", "Method"
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

        output$tableRevAll <- DT::renderDataTable(
            {
                shiny::req(selectedAssayList())
                x <- roundDbls(makeListTable(as.data.frame(selectedAssayList()[[2]])))
                names(x) <- c(
                    "Sequence", "GC content", "Tm", "Delta G"
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

        output$tablePr <- DT::renderDataTable(
            {
                shiny::req(length(selectedAssayList()) == 3)
                x <- roundDbls(removeListColumns(selectedAssayList()[[3]]))
                x$iupacSequence <- ifelse(x$plus, x$iupacSequence, NA)
                x$iupacSequenceRc <- ifelse(x$minus, x$iupacSequenceRc, NA)
                x <- x[, c(-1, -2)]
                names(x) <- c(
                    "Start", "End", "Length", "IUPAC sequence, plus",
                    "Iupac sequence, minus",
                    "Identity",
                    "Coverage", "Degeneracy", "GC, mean", "GC, range", "Tm, mean",
                    "Tm, range", "Delta G, mean", "Delta G, range", "Method"
                )
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

        output$tablePrAll <- DT::renderDataTable(
            {
                shiny::req(length(selectedAssayList()) == 3)
                x <- roundDbls(makeListTable(as.data.frame(selectedAssayList()[[3]])))
                names(x) <- c(
                    "Sequence, plus", "Sequence, minus", "GC content", "Tm", "Delta G"
                )
                if (!selectedAssayList()[[3]]$plus) {
                    x <- x[c("Sequence, minus", "GC content", "Tm", "Delta G")]
                    names(x) <- c(
                        "Sequence, minus", "GC content", "Tm", "Delta G"
                    )
                } else if (!selectedAssayList()[[3]]$minus) {
                    x <- x[c("Sequence, plus", "GC content", "Tm", "Delta G")]
                    names(x) <- c(
                        "Sequence, plus", "GC content", "Tm", "Delta G"
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

        output$tableFwdMatch <- DT::renderDataTable(
            {
                shiny::req(selectedAssayMatchList())
                x <- roundDbls(removeListColumns(selectedAssayMatchList()[[1]]))
                names(x) <- c(
                    "IUPAC sequence",
                    "Perfect match", "1 mismatch", "2 mismatches", "3 mismatches",
                    "4 or more mismatches"
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

        output$tableRevMatch <- DT::renderDataTable(
            {
                shiny::req(selectedAssayMatchList())
                x <- roundDbls(removeListColumns(selectedAssayMatchList()[[2]]))
                names(x) <- c(
                    "IUPAC sequence",
                    "Perfect match", "1 mismatch", "2 mismatches", "3 mismatches",
                    "4 or more mismatches"
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

        output$tablePrMatch <- DT::renderDataTable(
            {
                shiny::req(length(selectedAssayList()) == 3)
                x <- roundDbls(removeListColumns(selectedAssayMatchList()[[3]]))
                names(x) <- c(
                    "IUPAC sequence",
                    "Perfect match", "1 mismatch", "2 mismatches", "3 mismatches",
                    "4 or more mismatches"
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

        output$matchPlot <- shiny::renderPlot({
            shiny::req(is(selectedAssayMatch(), "RprimerMatchAssay"))
            plotData(selectedAssayMatch())
        })

        output$perfectMatchFwd <- shiny::renderText({
            shiny::req(!is.null(selectedAssayMatch()))
            c(
                "<b>Perfect match</b><br>",
                selectedAssayMatch()$idPerfectMatchFwd[[1]]
            )
        })


        output$oneMismatchFwd <- shiny::renderText({
            shiny::req(!is.null(selectedAssayMatch()))
            c(
                "<b>One mismatch</b><br>",
                selectedAssayMatch()$idOneMismatchFwd[[1]]
            )
        })

        output$twoMismatchesFwd <- shiny::renderText({
            shiny::req(!is.null(selectedAssayMatch()))
            c(
                "<b>Two mismatches</b><br>",
                selectedAssayMatch()$idTwoMismatchesFwd[[1]]
            )
        })

        output$threeMismatchesFwd <- shiny::renderText({
            shiny::req(!is.null(selectedAssayMatch()))
            c(
                "<b>Three mismatches</b><br>",
                selectedAssayMatch()$idThreeMismatchesFwd[[1]]
            )
        })

        output$fourOrMoreMismatchesFwd <- shiny::renderText({
            shiny::req(!is.null(selectedAssayMatch()))
            c(
                "<b>Four or more mismatches</b><br>",
                selectedAssayMatch()$idFourOrMoreMismatchesFwd[[1]]
            )
        })

        output$perfectMatchRev <- shiny::renderText({
            shiny::req(!is.null(selectedAssayMatch()))
            c(
                "<b>Perfect match</b><br>",
                selectedAssayMatch()$idPerfectMatchRev[[1]]
            )
        })

        output$oneMismatchRev <- shiny::renderText({
            shiny::req(!is.null(selectedAssayMatch()))
            c(
                "<b>One mismatch</b><br>",
                selectedAssayMatch()$idOneMismatchRev[[1]]
            )
        })

        output$twoMismatchesRev <- shiny::renderText({
            shiny::req(!is.null(selectedAssayMatch()))
            c(
                "<b>Two mismatches</b><br>",
                selectedAssayMatch()$idTwoMismatchesRev[[1]]
            )
        })

        output$threeMismatchesRev <- shiny::renderText({
            shiny::req(!is.null(selectedAssayMatch()))
            c(
                "<b>Three mismatches</b><br>",
                selectedAssayMatch()$idThreeMismatchesRev[[1]]
            )
        })

        output$fourOrMoreMismatchesRev <- shiny::renderText({
            shiny::req(!is.null(selectedAssayMatch()))
            c(
                "<b>Four or more mismatches</b><br>",
                selectedAssayMatch()$idFourOrMoreMismatchesRev[[1]]
            )
        })

        output$perfectMatchPr <- shiny::renderText({
            shiny::req(any(grepl("Pr", names(selectedAssayMatch()))))
            c(
                "<b>Perfect match</b><br>",
                selectedAssayMatch()$idPerfectMatchPr[[1]]
            )
        })

        output$oneMismatchPr <- shiny::renderText({
            shiny::req(any(grepl("Pr", names(selectedAssayMatch()))))
            c(
                "<b>One mismatch</b><br>",
                selectedAssayMatch()$idOneMismatchPr[[1]]
            )
        })

        output$twoMismatchesPr <- shiny::renderText({
            shiny::req(any(grepl("Pr", names(selectedAssayMatch()))))
            c(
                "<b>Two mismatches</b><br>",
                selectedAssayMatch()$idTwoMismatchesPr[[1]]
            )
        })

        output$threeMismatchesPr <- shiny::renderText({
            shiny::req(any(grepl("Pr", names(selectedAssayMatch()))))
            c(
                "<b>Three mismatches</b><br>",
                selectedAssayMatch()$idThreeMismatchesPr[[1]]
            )
        })

        output$fourOrMoreMismatchesPr <- shiny::renderText({
            shiny::req(any(grepl("Pr", names(selectedAssayMatch()))))
            c(
                "<b>Four or more mismatches</b><br>",
                selectedAssayMatch()$idFourOrMoreMismatchesPr[[1]]
            )
        })

        list(data = shiny::reactive(assayFilter()))
    })
}
