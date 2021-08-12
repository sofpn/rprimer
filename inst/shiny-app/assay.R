assayUI <- function(id) {
    ns <- NS(id)

    tagList(
        titlePanel("Step 4/5: Design assays"),
        br(),
        sidebarLayout(
            sidebarPanel(
                h4("Assay settings"),
                hr(),
                sliderInput(
                    ns("length"),
                    h5("Length"),
                    value = c(60, 120), min = 40, max = 5000,
                ),
                numericInput(
                    ns("tmDifferencePrimers"),
                    h5(
                        "Maximum melting temperature
                        difference between primers (Celcius degrees)"
                    ),
                    value = 10, min = 0, max = Inf,
                ),
                hr(),
                actionButton(ns("getAssays"), "Get assays")
            ),
            mainPanel(
                tabsetPanel(
                    id = ns("wizard"),
                    tabPanel(
                        title = "Plot",
                        br(),
                        shinycssloaders::withSpinner(
                            plotOutput(
                                ns("assayPlot"),
                                width = "100%",
                                height = 600
                            ),
                            color = "grey"
                        )
                    ),
                    tabPanel(
                        title = "Table",
                        br(),
                        uiOutput(ns("getData")),
                        br(),
                        h5("Select a row for more details"),
                        br(),
                        DT::dataTableOutput(ns("assayTable"))
                    ),
                    tabPanel(
                        title = "Selection",
                        br(),
                        width = 12,
                        tabsetPanel(
                            tabPanel(
                                title = "Overview",
                                br(),
                                uiOutput(ns("getSelectionData")),
                                br(),
                                h5("General assay information"),
                                hr(),
                                DT::dataTableOutput(ns("assayTableSelection")),
                                br(),
                                h5("Details"),
                                hr(),
                                br(),
                                uiOutput(ns("detailsTab"))
                            ),
                            tabPanel(
                                title = "Match details",
                                br(),
                                h5(
                                    "Proportion of matching sequences within the
                                                  intended target binding region in the input alignment"
                                ),
                                hr(),
                                br(),
                                column(
                                    width = 12, align = "center",
                                    shinycssloaders::withSpinner(
                                        plotOutput(
                                            ns("matchPlot"),
                                            width = "100%"
                                        ),
                                        color = "grey"
                                    )
                                ),
                                h5("Details"),
                                hr(),
                                br(),
                                uiOutput(ns("matchDetailsTab"))
                            )
                        )
                    )
                )
            )
        )
    )
}

assayServer <- function(id, alignment, consensus, oligos) {
    moduleServer(id, function(input, output, session) {

        assay <- eventReactive(input$getAssays, {
            req(is(oligos(), "RprimerOligo"))

            tryCatch(
                {
                    if (is.na(oligos()$length[[1]])) {
                        stop("No assays were found.", .call = FALSE)
                    }

                    assays(oligos(),
                        length = input$length,
                        tmDifferencePrimers = as.numeric(input$tmDifferencePrimers)
                    )

                },
                error = function(cond) {
                    showNotification(
                        "No assays were found.\n
                        Try to adjust design settings.",
                        type = "error",
                        duration = NULL
                    )
                },
                silent = TRUE
            )
        })

        output$assayPlot <- renderPlot({
            req(is(assay(), "RprimerAssay"))
            plotData(assay())
        })

        observeEvent(input$assayTable_rows_selected, {
            Sys.sleep(1)
            updateTabsetPanel(session, "wizard", selected = "Selection")
        })

        selectedAssay <- reactive({
            req(is(assay(), "RprimerAssay"))
            if (!is.null(input$assayTable_rows_selected)) {
                assay()[input$assayTable_rows_selected, ]
            } else {
                NULL
            }
        })

        selectedAssayList <- reactive({
            req(selectedAssay())
            splitAssayToList(selectedAssay())
        })

        selectedAssayMatch <- reactive({
            req(selectedAssay())
            req(alignment())
            if (is.na(selectedAssay()$length[[1]])) {
                NULL
            } else {
                checkMatch(selectedAssay(), alignment())
            }
        })

        selectedAssayMatchList <- reactive({
            req(selectedAssayMatch())
            splitAssayToList(selectedAssayMatch())
        })

        output$downloadTable <- downloadHandler(
            filename <- function() {
                paste0("assays-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                write.table(
                    as.data.frame(assaySelection()), file,
                    quote = FALSE, sep = "\t",
                    row.names = FALSE
                )
            }
        )

        output$getData <- renderUI({
            req(is(assay(), "RprimerAssay"))
            ns <- session$ns
            fluidRow(
                downloadLink(
                    ns("downloadTable"), "Download table as .txt"
                ),
                br(),
                downloadLink(
                    ns("downloadFasta"),
                    "Download assays in fasta-format"
                )
            )
        })

        output$downloadFasta <- downloadHandler(
            filename <- function() {
                paste0("assays-fasta", Sys.Date(), ".txt")
            },
            content <- function(file) {
                x <- assayection()
                x <- as(assaySelection(), "DNAStringSet")
                Biostrings::writeXStringSet(x, file)
            }
        )


        output$assayTable <- DT::renderDataTable(
            {
                req(is(assay(), "RprimerAssay"))
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

        output$detailsTab <- renderUI({
            ns <- session$ns
            tabs <- list(
                tabPanel(
                    title = "Forward",
                    br(),
                    DT::dataTableOutput(ns("tableFwd")),
                    br(),
                    br(),
                    h5("All sequence variants"),
                    hr(),
                    br(),
                    DT::dataTableOutput(ns("tableFwdAll")),
                    br(),
                    br(),
                    h5("Nucleotide distribution in target alignment"),
                    hr(),
                    br(),
                    column(
                        width = 12, align = "center",
                        shinycssloaders::withSpinner(
                            plotOutput(
                                ns("ntPlotFwd"),
                                width = "75%"
                            ),
                            color = "grey"
                        )
                    ),
                    br(),
                    br()
                ),
                tabPanel(
                    title = "Reverse",
                    br(),
                    DT::dataTableOutput(ns("tableRev")),
                    br(),
                    br(),
                    h5("All sequence variants"),
                    hr(),
                    br(),
                    DT::dataTableOutput(ns("tableRevAll")),
                    br(),
                    br(),
                    h5("Nucleotide distribution in target alignment"),
                    hr(),
                    br(),
                    column(
                        width = 12, align = "center",
                        shinycssloaders::withSpinner(
                            plotOutput(
                                ns("ntPlotRev"),
                                width = "75%"
                            ),
                            color = "grey"
                        )
                    )
                ),
                tabPanel(
                    title = "Probe",
                    br(),
                    DT::dataTableOutput(ns("tablePr")),
                    br(),
                    br(),
                    h5("All sequence variants"),
                    hr(),
                    br(),
                    DT::dataTableOutput(ns("tablePrAll")),
                    br(),
                    br(),
                    h5("Nucleotide distribution in target alignment"),
                    hr(),
                    br(),
                    column(
                        width = 12, align = "center",
                        shinycssloaders::withSpinner(
                            plotOutput(
                                ns("ntPlotPr"),
                                width = "75%"
                            ),
                            color = "grey"
                        )
                    ),
                    br(),
                    br()
                )
            )

            if (length(selectedAssayList()) == 3) {
                do.call(tabsetPanel, tabs)
            } else {
                do.call(tabsetPanel, tabs[1:2])
            }
        })

        output$matchDetailsTab <- renderUI({
            ns <- session$ns

            tabs <- list(
                tabPanel(
                    title = "Forward",
                    br(),
                    DT::dataTableOutput(
                        ns("tableFwdMatch")
                    ),
                    br(),
                    h5("Target sequence names"),
                    hr(),
                    htmlOutput(ns("perfectMatchFwd")),
                    br(),
                    br(),
                    htmlOutput(ns("oneMismatchFwd")),
                    br(),
                    br(),
                    htmlOutput(ns("twoMismatchesFwd")),
                    br(),
                    br(),
                    htmlOutput(ns("threeMismatchesFwd")),
                    br(),
                    br(),
                    htmlOutput(ns("fourOrMoreMismatchesFwd")),
                    br(),
                    br()
                ),
                tabPanel(
                    title = "Reverse",
                    br(),
                    DT::dataTableOutput(
                        ns("tableRevMatch")
                    ),
                    br(),
                    h5("Target sequence names"),
                    hr(),
                    htmlOutput(ns("perfectMatchRev")),
                    br(),
                    br(),
                    htmlOutput(ns("oneMismatchRev")),
                    br(),
                    br(),
                    htmlOutput(ns("twoMismatchesRev")),
                    br(),
                    br(),
                    htmlOutput(ns("threeMismatchesRev")),
                    br(),
                    br(),
                    htmlOutput(ns("fourOrMoreMismatchesRev")),
                    br(),
                    br()
                ),
                tabPanel(
                    title = "Probe",
                    br(),
                    DT::dataTableOutput(
                        ns("tablePrMatch")
                    ),
                    br(),
                    h5("Target sequence names"),
                    hr(),
                    htmlOutput(ns("perfectMatchPr")),
                    br(),
                    br(),
                    htmlOutput(ns("oneMismatchPr")),
                    br(),
                    br(),
                    htmlOutput(ns("twoMismatchesPr")),
                    br(),
                    br(),
                    htmlOutput(ns("threeMismatchesPr")),
                    br(),
                    br(),
                    htmlOutput(ns("fourOrMoreMismatchesPr")),
                    br(),
                    br()
                )
            )

            if (length(selectedAssayList()) == 3) {
                do.call(tabsetPanel, tabs)
            } else {
                do.call(tabsetPanel, tabs[1:2])
            }
        })

        output$assayTableSelection <- DT::renderDataTable(
            {
                req(selectedAssay())
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

        output$ntPlotFwd <- renderPlot({
            req(selectedAssay())
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

        output$ntPlotRev <- renderPlot({
            req(selectedAssay())
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

        output$ntPlotPr <- renderPlot({
            req(selectedAssay())
            req(grepl("Pr", names(selectedAssay())))
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

        output$getSelectionData <- renderUI({
            req(selectedAssay())
            ns <- session$ns
            fluidRow(
                downloadLink(
                    ns("downloadSelectionTable"),
                    "Download assay information as .txt"
                ),
                br(),
                downloadLink(
                    ns("downloadSelectionFasta"),
                    "Download assay in fasta-format"
                )
            )
        })

        output$downloadSelectionTable <- downloadHandler(
            filename <- function() {
                paste0("assay-selection-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                write.table(
                    as.data.frame(selectedAssay()), file,
                    quote = FALSE, sep = "\t",
                    row.names = FALSE
                )
            }
        )

        output$downloadSelectionFasta <- downloadHandler(
            filename <- function() {
                paste0("assay-selection-fasta-", Sys.Date(), ".txt")
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
                req(selectedAssayList())
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
                req(selectedAssayList())
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
                req(selectedAssayList())
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
                req(length(selectedAssayList()) == 3)
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
                req(length(selectedAssayList()) == 3)
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
                req(selectedAssayMatchList())
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
                req(selectedAssayMatchList())
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
                req(length(selectedAssayList()) == 3)
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

        output$matchPlot <- renderPlot({
            req(is(selectedAssayMatch(), "RprimerMatchAssay"))
            plotData(selectedAssayMatch())
        })

        output$perfectMatchFwd <- renderText({
            req(!is.null(selectedAssayMatch()))
            c(
                "<b>Perfect match</b><br>",
                selectedAssayMatch()$idPerfectMatchFwd[[1]]
            )
        })


        output$oneMismatchFwd <- renderText({
            req(!is.null(selectedAssayMatch()))
            c(
                "<b>One mismatch</b><br>",
                selectedAssayMatch()$idOneMismatchFwd[[1]]
            )
        })

        output$twoMismatchesFwd <- renderText({
            req(!is.null(selectedAssayMatch()))
            c(
                "<b>Two mismatches</b><br>",
                selectedAssayMatch()$idTwoMismatchesFwd[[1]]
            )
        })

        output$threeMismatchesFwd <- renderText({
            req(!is.null(selectedAssayMatch()))
            c(
                "<b>Three mismatches</b><br>",
                selectedAssayMatch()$idThreeMismatchesFwd[[1]]
            )
        })

        output$fourOrMoreMismatchesFwd <- renderText({
            req(!is.null(selectedAssayMatch()))
            c(
                "<b>Four or more mismatches</b><br>",
                selectedAssayMatch()$idFourOrMoreMismatchesFwd[[1]]
            )
        })

        output$perfectMatchRev <- renderText({
            req(!is.null(selectedAssayMatch()))
            c(
                "<b>Perfect match</b><br>",
                selectedAssayMatch()$idPerfectMatchRev[[1]]
            )
        })

        output$oneMismatchRev <- renderText({
            req(!is.null(selectedAssayMatch()))
            c(
                "<b>One mismatch</b><br>",
                selectedAssayMatch()$idOneMismatchRev[[1]]
            )
        })

        output$twoMismatchesRev <- renderText({
            req(!is.null(selectedAssayMatch()))
            c(
                "<b>Two mismatches</b><br>",
                selectedAssayMatch()$idTwoMismatchesRev[[1]]
            )
        })

        output$threeMismatchesRev <- renderText({
            req(!is.null(selectedAssayMatch()))
            c(
                "<b>Three mismatches</b><br>",
                selectedAssayMatch()$idThreeMismatchesRev[[1]]
            )
        })

        output$fourOrMoreMismatchesRev <- renderText({
            req(!is.null(selectedAssayMatch()))
            c(
                "<b>Four or more mismatches</b><br>",
                selectedAssayMatch()$idFourOrMoreMismatchesRev[[1]]
            )
        })

        output$perfectMatchPr <- renderText({
            req(any(grepl("Pr", names(selectedAssayMatch()))))
            c(
                "<b>Perfect match</b><br>",
                selectedAssayMatch()$idPerfectMatchPr[[1]]
            )
        })

        output$oneMismatchPr <- renderText({
            req(any(grepl("Pr", names(selectedAssayMatch()))))
            c(
                "<b>One mismatch</b><br>",
                selectedAssayMatch()$idOneMismatchPr[[1]]
            )
        })

        output$twoMismatchesPr <- renderText({
            req(any(grepl("Pr", names(selectedAssayMatch()))))
            c(
                "<b>Two mismatches</b><br>",
                selectedAssayMatch()$idTwoMismatchesPr[[1]]
            )
        })

        output$threeMismatchesPr <- renderText({
            req(any(grepl("Pr", names(selectedAssayMatch()))))
            c(
                "<b>Three mismatches</b><br>",
                selectedAssayMatch()$idThreeMismatchesPr[[1]]
            )
        })

        output$fourOrMoreMismatchesPr <- renderText({
            req(any(grepl("Pr", names(selectedAssayMatch()))))
            c(
                "<b>Four or more mismatches</b><br>",
                selectedAssayMatch()$idFourOrMoreMismatchesPr[[1]]
            )
        })

        list(data = reactive(assay()))
    })
}
