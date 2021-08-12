oligoFilterUI <- function(id) {
    ns <- NS(id)

    tagList(
        titlePanel("Step 3/5: Filter oligos"),
        br(),
        sidebarLayout(
            sidebarPanel(
                uiOutput(ns("fwd")),
                uiOutput(ns("rev")),
                uiOutput(ns("pr"))
            ),
            mainPanel(
                tabsetPanel(
                    id = ns("wizard"),
                    tabPanel(
                        title = "Plot",
                        br(),
                        shinycssloaders::withSpinner(
                            plotOutput(
                                ns("oligoPlot"),
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
                        DT::dataTableOutput(ns("oligoTable"))
                    ),
                    tabPanel(
                        title = "Selection",
                        br(),
                        tabsetPanel(
                            tabPanel(
                                title = "Overview",
                                br(),
                                uiOutput(ns("getSelectionData")),
                                br(),
                                h5("Oligo information"),
                                hr(),
                                DT::dataTableOutput(ns("oligoTableSelection")),
                                br(),
                                h5("All sequence variants"),
                                hr(),
                                DT::dataTableOutput(ns("oligoTableSelectionAll")),
                                br(),
                                h5("Nucleotide distribution in target alignment (5'-3')"),
                                hr(),
                                br(),
                                column(
                                    width = 12, align = "center",
                                    shinycssloaders::withSpinner(
                                        plotOutput(
                                            ns("ntPlot"),
                                            width = "75%"
                                        ),
                                        color = "grey"
                                    )
                                )
                            ),
                            tabPanel(
                                title = "Match details",
                                br(),
                                h5(
                                    "Proportion of matching sequences
                                        within the intended target binding
                                        region in the input alignmnent"
                                ),
                                hr(),
                                DT::dataTableOutput(
                                    ns("oligoTableSelectionMatch")
                                ),
                                br(),
                                br(),
                                column(
                                    width = 12, align = "center",
                                    shinycssloaders::withSpinner(
                                        plotOutput(
                                            ns("matchPlot"),
                                            width = "75%"
                                        ),
                                        color = "grey"
                                    )
                                ),
                                br(),
                                br(),
                                h5("Target sequence names"),
                                hr(),
                                htmlOutput(ns("perfectMatch")),
                                br(),
                                br(),
                                htmlOutput(ns("oneMismatch")),
                                br(),
                                br(),
                                htmlOutput(ns("twoMismatches")),
                                br(),
                                br(),
                                htmlOutput(ns("threeMismatches")),
                                br(),
                                br(),
                                htmlOutput(ns("fourOrMoreMismatches")),
                                br(),
                                br()
                            )
                        )
                    )
                )
            )
        )
    )
}

oligoFilterServer <- function(id, alignment, consensus, oligos) {
    moduleServer(id, function(input, output, session) {
        output$fwd <- renderUI({
            req(any(oligos()$type == "primer" & oligos()$fwd))

            ns <- session$ns

            list(
                h5("Forward"),
                hr(),
                numericInput(ns("fwdRegionFrom"), h5("From"),
                    min = oligos()$roiStart[[1]],
                    max = oligos()$roiEnd[[1]],
                    value = oligos()$roiStart[[1]]
                ),
                numericInput(ns("fwdRegionTo"), h5("To"),
                    min = oligos()$roiStart[[1]],
                    max = oligos()$roiEnd[[1]],
                    value = oligos()$roiEnd[[1]]
                ),
                sliderInput(ns("minOligoIdentityFwd"),
                    round = -4, step = 0.0001,
                    h5("Mminimum identity"),
                    min = round(min(
                        oligos()$identity[
                            oligos()$fwd &
                                oligos()$type == "primer"
                        ],
                        na.rm = TRUE
                    ), 4),
                    max = round(max(
                        oligos()$identity[
                            oligos()$fwd &
                                oligos()$type == "primer"
                        ],
                        na.rm = TRUE
                    ), 4),
                    value = round(min(
                        oligos()$identity[
                            oligos()$fwd &
                                oligos()$type == "primer"
                        ],
                        na.rm = TRUE
                    ), 4)
                ),
                sliderInput(ns("minOligoCoverageFwd"),
                    round = -4, step = 0.0001,
                    h5("Minimum coverage"),
                    min = round(min(
                        oligos()$coverage[
                            oligos()$fwd &
                                oligos()$type == "primer"
                        ],
                        na.rm = TRUE
                    ), 4),
                    max = round(max(
                        oligos()$coverage[
                            oligos()$fwd &
                                oligos()$type == "primer"
                        ],
                        na.rm = TRUE
                    ), 4),
                    value = round(min(
                        oligos()$coverage[
                            oligos()$fwd &
                                oligos()$type == "primer"
                        ],
                        na.rm = TRUE
                    ), 4)
                )
            )
        })

        output$rev <- renderUI({
            req(any(oligos()$type == "primer" & oligos()$rev))

            ns <- session$ns

            list(
                h5("Reverse"),
                hr(),
                numericInput(ns("revRegionFrom"), h5("From"),
                    min = oligos()$roiStart[[1]],
                    max = oligos()$roiEnd[[1]],
                    value = oligos()$roiStart[[1]]
                ),
                numericInput(ns("revRegionTo"), h5("To"),
                    min = oligos()$roiStart[[1]],
                    max = oligos()$roiEnd[[1]],
                    value = oligos()$roiEnd[[1]]
                ),
                sliderInput(ns("minOligoIdentityRev"),
                    round = -4, step = 0.0001,
                    h5("Minimum identity"),
                    min = round(min(
                        oligos()$identity[
                            oligos()$rev &
                                oligos()$type == "primer"
                        ],
                        na.rm = TRUE
                    ), 4),
                    max = round(max(
                        oligos()$identity[
                            oligos()$rev &
                                oligos()$type == "primer"
                        ],
                        na.rm = TRUE
                    ), 4),
                    value = round(min(
                        oligos()$identity[
                            oligos()$rev &
                                oligos()$type == "primer"
                        ],
                        na.rm = TRUE
                    ), 4)
                ),
                sliderInput(ns("minOligoCoverageRev"),
                    round = -4, step = 0.0001,
                    h5("Minimum coverage"),
                    min = round(min(
                        oligos()$coverage[
                            oligos()$rev &
                                oligos()$type == "primer"
                        ],
                        na.rm = TRUE
                    ), 4),
                    max = round(max(
                        oligos()$coverage[
                            oligos()$rev &
                                oligos()$type == "primer"
                        ],
                        na.rm = TRUE
                    ), 4),
                    value = round(min(
                        oligos()$coverage[
                            oligos()$rev &
                                oligos()$type == "primer"
                        ],
                        na.rm = TRUE
                    ), 4)
                )
            )
        })


        output$pr <- renderUI({
            req(any(oligos()$type == "probe"))

            ns <- session$ns

            list(
                h5("Probe"),
                hr(),
                numericInput(ns("prRegionFrom"), h5("From"),
                    min = oligos()$roiStart[[1]],
                    max = oligos()$roiEnd[[1]],
                    value = oligos()$roiStart[[1]]
                ),
                numericInput(ns("prRegionTo"), h5("To"),
                    min = oligos()$roiStart[[1]],
                    max = oligos()$roiEnd[[1]],
                    value = oligos()$roiEnd[[1]]
                ),
                sliderInput(ns("minOligoIdentityPr"),
                    round = -4, step = 0.0001,
                    h5("Minimum identity"),
                    min = round(min(
                        oligos()$identity[
                            oligos()$type == "probe"
                        ],
                        na.rm = TRUE
                    ), 4),
                    max = round(max(
                        oligos()$identity[
                            oligos()$type == "probe"
                        ],
                        na.rm = TRUE
                    ), 4),
                    value = round(min(
                        oligos()$identity[
                            oligos()$type == "probe"
                        ],
                        na.rm = TRUE
                    ), 4)
                ),
                sliderInput(ns("minOligoCoveragePr"),
                    round = -4, step = 0.0001,
                    h5("Minimum coverage"),
                    min = round(min(
                        oligos()$coverage[
                            oligos()$type == "probe"
                        ],
                        na.rm = TRUE
                    ), 4),
                    max = round(max(
                        oligos()$coverage[
                            oligos()$type == "probe"
                        ],
                        na.rm = TRUE
                    ), 4),
                    value = round(min(
                        oligos()$coverage[
                            oligos()$type == "probe"
                        ],
                        na.rm = TRUE
                    ), 4)
                )
            )
        })


        ols <- reactive({
            req(oligos())
            filterOligos(oligos(),
                fwdFrom = input$fwdRegionFrom,
                fwdTo = input$fwdRegionTo,
                revFrom = input$revRegionFrom,
                revTo = input$revRegionTo,
                prFrom = input$prRegionFrom,
                prTo = input$prRegionTo,
                fwdIdentity = input$minOligoIdentityFwd,
                revIdentity = input$minOligoIdentityRev,
                prIdentity = input$minOligoIdentityPr,
                fwdCoverage = input$minOligoCoverageFwd,
                revCoverage = input$minOligoCoverageRev,
                prCoverage = input$minOligoCoveragePr
            )
        })

        selectedOligo <- reactive({
            req(ols())
            if (!is.null(input$oligoTable_rows_selected)) {
                ols()[input$oligoTable_rows_selected, ]
            } else {
                NULL
            }
        })

        observeEvent(input$oligoTable_rows_selected, {
            Sys.sleep(1)
            updateTabsetPanel(session, "wizard", selected = "Selection")
        })

        selectedOligoMatch <- reactive({
            req(is(selectedOligo(), "RprimerOligo"))
            req(is(alignment(), "DNAMultipleAlignment"))
            checkMatch(selectedOligo(), alignment())
        })

        output$oligoPlot <- renderPlot({
            req(is(ols(), "RprimerOligo"))
            plotData(ols())
        })

        output$getData <- renderUI({
            req(!is.na(ols()$length[[1]]))
            ns <- session$ns
            fluidRow(
                downloadLink(
                    ns("downloadTable"), "Download table as .txt"
                ),
                br(),
                downloadLink(
                    ns("downloadFasta"),
                    "Download oligo sequences in fasta-format"
                )
            )
        })

        output$downloadTable <- downloadHandler(
            filename <- function() {
                paste0("oligos-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                write.table(
                    as.data.frame(ols()), file,
                    quote = FALSE, sep = "\t",
                    row.names = FALSE
                )
            }
        )

        output$downloadFasta <- downloadHandler(
            filename <- function() {
                paste0("oligos-fasta", Sys.Date(), ".txt")
            },
            content <- function(file) {
                x <- ols()
                x <- as(ols(), "DNAStringSet")
                Biostrings::writeXStringSet(x, file)
            }
        )

        output$oligoTable <- DT::renderDataTable(
            {
                req(is(ols(), "RprimerOligo"))
                if (is.na(ols()$length[[1]])) {
                    NULL
                } else {
                    x <- roundDbls(removeListColumns(as.data.frame(ols())))
                    names(x) <- c(
                        "Type", "Forward", "Reverse", "Start", "End", "Length",
                        "IUPAC sequence", "IUPAC sequence, RC", "Identity",
                        "Coverage", "Degeneracy", "GC content, mean",
                        "GC content, range", "Tm, mean", "Tm, range", "Delta G, mean",
                        "Delta G, range", "Design method", "Score", "ROI, start",
                        "ROI, end"
                    )
                    x
                }
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

        output$getSelectionData <- renderUI({
            req(is(selectedOligo(), "RprimerOligo"))
            ns <- session$ns
            fluidRow(
                downloadLink(
                    ns("downloadSelectionTable"),
                    "Download oligo information as .txt"
                ),
                br(),
                downloadLink(
                    ns("downloadSelectionFasta"),
                    "Download oligo sequence in fasta-format"
                )
            )
        })

        output$downloadSelectionTable <- downloadHandler(
            filename <- function() {
                paste0("oligo-selection-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                write.table(
                    as.data.frame(selectedOligo()), file,
                    quote = FALSE, sep = "\t",
                    row.names = FALSE
                )
            }
        )

        output$downloadSelectionFasta <- downloadHandler(
            filename <- function() {
                paste0("oligo-selection-fasta-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                x <- selectedOligo()
                x <- as(selectedOligo(), "DNAStringSet")
                Biostrings::writeXStringSet(x, file)
            }
        )

        output$oligoTableSelection <- DT::renderDataTable(
            {
                req(is(selectedOligo(), "RprimerOligo"))
                x <- roundDbls(removeListColumns(as.data.frame(selectedOligo())))
                names(x) <- c(
                    "Type", "Forward", "Reverse", "Start", "End", "Length",
                    "IUPAC sequence", "IUPAC sequence, RC", "Identity",
                    "Coverage", "Degeneracy", "GC content, mean",
                    "GC content, range", "Tm, mean", "Tm, range", "Delta G, mean",
                    "Delta G, range", "Design method", "Score", "ROI, start",
                    "ROI, end"
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


        output$oligoTableSelectionAll <- DT::renderDataTable(
            {
                req(is(selectedOligo(), "RprimerOligo"))
                x <- roundDbls(makeListTable(as.data.frame(selectedOligo())))
                names(x) <- c(
                    "Sequence", "Sequence, RC", "GC content", "Tm", "Delta G"
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


        output$ntPlot <- renderPlot({
            req(is(selectedOligo(), "RprimerOligo"))
            from <- selectedOligo()$start
            to <- selectedOligo()$end
            plotData(consensus()[
                consensus()$position >= from & consensus()$position <= to,
            ], type = "nucleotide")
        })

        output$oligoTableSelectionMatch <- DT::renderDataTable(
            {
                req(is(selectedOligoMatch(), "RprimerMatchOligo"))
                x <- roundDbls(removeListColumns(as.data.frame(selectedOligoMatch())))
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

        output$matchPlot <- renderPlot({
            req(is(selectedOligoMatch(), "RprimerMatchOligo"))
            plotData(selectedOligoMatch())
        })

        output$perfectMatch <- renderText({
            req(is(selectedOligoMatch(), "RprimerMatchOligo"))
            c(
                "<b>Perfect match</b><br>",
                selectedOligoMatch()$idPerfectMatch[[1]]
            )
        })

        output$oneMismatch <- renderText({
            req(is(selectedOligoMatch(), "RprimerMatchOligo"))
            c(
                "<b>One mismatch</b><br>",
                selectedOligoMatch()$idOneMismatch[[1]]
            )
        })

        output$twoMismatches <- renderText({
            req(is(selectedOligoMatch(), "RprimerMatchOligo"))
            c(
                "<b>Two mismatches</b><br>",
                selectedOligoMatch()$idTwoMismatches[[1]]
            )
        })

        output$threeMismatches <- renderText({
            req(is(selectedOligoMatch(), "RprimerMatchOligo"))
            c(
                "<b>Three mismatches</b><br>",
                selectedOligoMatch()$idThreeMismatches[[1]]
            )
        })

        output$fourOrMoreMismatches <- renderText({
            req(is(selectedOligoMatch(), "RprimerMatchOligo"))
            c(
                "<b>Four or more mismatches</b><br>",
                selectedOligoMatch()$idFourOrMoreMismatches[[1]]
            )
        })

        list(data = reactive(ols()))
    })
}
