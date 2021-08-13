oligoFilterUI <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::titlePanel("Step 3/5: Filter oligos"),
        shiny::br(),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::uiOutput(ns("fwd")),
                shiny::uiOutput(ns("rev")),
                shiny::uiOutput(ns("pr"))
            ),
            shiny::mainPanel(
                shiny::tabsetPanel(
                    id = ns("wizard"),
                    shiny::tabPanel(
                        title = "Plot",
                        shiny::br(),
                        shinycssloaders::withSpinner(
                            shiny::plotOutput(
                                ns("oligoPlot"),
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
                        DT::dataTableOutput(ns("oligoTable"))
                    ),
                    shiny::tabPanel(
                        title = "Selection",
                        shiny::br(),
                        shiny::tabsetPanel(
                            shiny::tabPanel(
                                title = "Overview",
                                shiny::br(),
                                shiny::uiOutput(ns("getSelectionData")),
                                shiny::br(),
                                shiny::h5("Oligo information"),
                                shiny::hr(),
                                DT::dataTableOutput(ns("oligoTableSelection")),
                                shiny::br(),
                                shiny::h5("All sequence variants"),
                                shiny::hr(),
                                DT::dataTableOutput(ns("oligoTableSelectionAll")),
                                shiny::br(),
                                shiny::h5("Nucleotide distribution in target alignment (5'-3')"),
                                shiny::hr(),
                                shiny::br(),
                                shiny::column(
                                    width = 12, align = "center",
                                    shinycssloaders::withSpinner(
                                        shiny::plotOutput(
                                            ns("ntPlot"),
                                            width = "75%"
                                        ),
                                        color = "grey"
                                    )
                                )
                            ),
                            shiny::tabPanel(
                                title = "Match details",
                                shiny::br(),
                                shiny::h5(
                                    "Proportion of matching sequences
                                        within the intended target binding
                                        region in the input alignmnent"
                                ),
                                shiny::hr(),
                                DT::dataTableOutput(
                                    ns("oligoTableSelectionMatch")
                                ),
                                shiny::br(),
                                shiny::br(),
                                shiny::column(
                                    width = 12, align = "center",
                                    shinycssloaders::withSpinner(
                                        shiny::plotOutput(
                                            ns("matchPlot"),
                                            width = "75%"
                                        ),
                                        color = "grey"
                                    )
                                ),
                                shiny::br(),
                                shiny::br(),
                                shiny::h5("Target sequence names"),
                                shiny::hr(),
                                shiny::htmlOutput(ns("perfectMatch")),
                                shiny::br(),
                                shiny::br(),
                                shiny::htmlOutput(ns("oneMismatch")),
                                shiny::br(),
                                shiny::br(),
                                shiny::htmlOutput(ns("twoMismatches")),
                                shiny::br(),
                                shiny::br(),
                                shiny::htmlOutput(ns("threeMismatches")),
                                shiny::br(),
                                shiny::br(),
                                shiny::htmlOutput(ns("fourOrMoreMismatches")),
                                shiny::br(),
                                shiny::br()
                            )
                        )
                    )
                )
            )
        )
    )
}

oligoFilterServer <- function(id, alignment, consensus, oligos) {
    shiny::moduleServer(id, function(input, output, session) {
        output$fwd <- shiny::renderUI({
            shiny::req(any(oligos()$type == "primer" & oligos()$fwd))

            ns <- session$ns

            list(
                shiny::h5("Forward"),
                shiny::hr(),
                shiny::numericInput(ns("fwdRegionFrom"), shiny::h5("From"),
                    min = oligos()$roiStart[[1]],
                    max = oligos()$roiEnd[[1]],
                    value = oligos()$roiStart[[1]]
                ),
                shiny::numericInput(ns("fwdRegionTo"), shiny::h5("To"),
                    min = oligos()$roiStart[[1]],
                    max = oligos()$roiEnd[[1]],
                    value = oligos()$roiEnd[[1]]
                ),
                shiny::sliderInput(ns("minOligoIdentityFwd"),
                    round = -4, step = 0.0001,
                    shiny::h5("Mminimum identity"),
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
                shiny::sliderInput(ns("minOligoCoverageFwd"),
                    round = -4, step = 0.0001,
                    shiny::h5("Minimum coverage"),
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

        output$rev <- shiny::renderUI({
            shiny::req(any(oligos()$type == "primer" & oligos()$rev))

            ns <- session$ns

            list(
                shiny::h5("Reverse"),
                shiny::hr(),
                shiny::numericInput(ns("revRegionFrom"), shiny::h5("From"),
                    min = oligos()$roiStart[[1]],
                    max = oligos()$roiEnd[[1]],
                    value = oligos()$roiStart[[1]]
                ),
                shiny::numericInput(ns("revRegionTo"), shiny::h5("To"),
                    min = oligos()$roiStart[[1]],
                    max = oligos()$roiEnd[[1]],
                    value = oligos()$roiEnd[[1]]
                ),
                shiny::sliderInput(ns("minOligoIdentityRev"),
                    round = -4, step = 0.0001,
                    shiny::h5("Minimum identity"),
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
                shiny::sliderInput(ns("minOligoCoverageRev"),
                    round = -4, step = 0.0001,
                    shiny::h5("Minimum coverage"),
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


        output$pr <- shiny::renderUI({
            shiny::req(any(oligos()$type == "probe"))

            ns <- session$ns

            list(
                shiny::h5("Probe"),
                shiny::hr(),
                shiny::numericInput(ns("prRegionFrom"), shiny::h5("From"),
                    min = oligos()$roiStart[[1]],
                    max = oligos()$roiEnd[[1]],
                    value = oligos()$roiStart[[1]]
                ),
                shiny::numericInput(ns("prRegionTo"), shiny::h5("To"),
                    min = oligos()$roiStart[[1]],
                    max = oligos()$roiEnd[[1]],
                    value = oligos()$roiEnd[[1]]
                ),
                shiny::sliderInput(ns("minOligoIdentityPr"),
                    round = -4, step = 0.0001,
                    shiny::h5("Minimum identity"),
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
                shiny::sliderInput(ns("minOligoCoveragePr"),
                    round = -4, step = 0.0001,
                    shiny::h5("Minimum coverage"),
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


        ols <- shiny::reactive({
            shiny::req(oligos())
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

        selectedOligo <- shiny::reactive({
            shiny::req(ols())
            if (!is.null(input$oligoTable_rows_selected)) {
                ols()[input$oligoTable_rows_selected, ]
            } else {
                NULL
            }
        })

        shiny::observeEvent(input$oligoTable_rows_selected, {
            Sys.sleep(1)
            shiny::updateTabsetPanel(session, "wizard", selected = "Selection")
        })

        selectedOligoMatch <- shiny::reactive({
            shiny::req(is(selectedOligo(), "RprimerOligo"))
            shiny::req(is(alignment(), "DNAMultipleAlignment"))
            checkMatch(selectedOligo(), alignment())
        })

        output$oligoPlot <- shiny::renderPlot({
            shiny::req(is(ols(), "RprimerOligo"))
            plotData(ols())
        })

        output$getData <- shiny::renderUI({
            shiny::req(!is.na(ols()$length[[1]]))
            ns <- session$ns
            list(
                shiny::downloadLink(
                    ns("downloadTable"), "Download table as .txt"
                ),
                shiny::br(),
                shiny::downloadLink(
                    ns("downloadFasta"),
                    "Download oligo sequences in fasta-format"
                )
            )
        })

        output$downloadTable <- shiny::downloadHandler(
            filename <- function() {
                paste0("oligos-filtered", Sys.Date(), ".txt")
            },
            content <- function(file) {
                write.table(
                    as.data.frame(ols()), file,
                    quote = FALSE, sep = "\t",
                    row.names = FALSE
                )
            }
        )

        output$downloadFasta <- shiny::downloadHandler(
            filename <- function() {
                paste0("oligos-filtered-fasta", Sys.Date(), ".txt")
            },
            content <- function(file) {
                x <- ols()
                x <- as(ols(), "DNAStringSet")
                Biostrings::writeXStringSet(x, file)
            }
        )

        output$oligoTable <- DT::renderDataTable(
            {
                shiny::req(is(ols(), "RprimerOligo"))
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

        output$getSelectionData <- shiny::renderUI({
            shiny::req(is(selectedOligo(), "RprimerOligo"))
            ns <- session$ns
            list(
                shiny::downloadLink(
                    ns("downloadSelectionTable"),
                    "Download oligo information as .txt"
                ),
                shiny::br(),
                shiny::downloadLink(
                    ns("downloadSelectionFasta"),
                    "Download oligo sequence in fasta-format"
                )
            )
        })

        output$downloadSelectionTable <- shiny::downloadHandler(
            filename <- function() {
                paste0("oligo-filtered-selection-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                write.table(
                    as.data.frame(selectedOligo()), file,
                    quote = FALSE, sep = "\t",
                    row.names = FALSE
                )
            }
        )

        output$downloadSelectionFasta <- shiny::downloadHandler(
            filename <- function() {
                paste0("oligo-filtered-selection-fasta-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                x <- selectedOligo()
                x <- as(selectedOligo(), "DNAStringSet")
                Biostrings::writeXStringSet(x, file)
            }
        )

        output$oligoTableSelection <- DT::renderDataTable(
            {
                shiny::req(is(selectedOligo(), "RprimerOligo"))
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
                shiny::req(is(selectedOligo(), "RprimerOligo"))
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


        output$ntPlot <- shiny::renderPlot({
            shiny::req(is(selectedOligo(), "RprimerOligo"))
            from <- selectedOligo()$start
            to <- selectedOligo()$end
            plotData(consensus()[
                consensus()$position >= from & consensus()$position <= to,
            ], type = "nucleotide")
        })

        output$oligoTableSelectionMatch <- DT::renderDataTable(
            {
                shiny::req(is(selectedOligoMatch(), "RprimerMatchOligo"))
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

        output$matchPlot <- shiny::renderPlot({
            shiny::req(is(selectedOligoMatch(), "RprimerMatchOligo"))
            plotData(selectedOligoMatch())
        })

        output$perfectMatch <- shiny::renderText({
            shiny::req(is(selectedOligoMatch(), "RprimerMatchOligo"))
            c(
                "<b>Perfect match</b><br>",
                selectedOligoMatch()$idPerfectMatch[[1]]
            )
        })

        output$oneMismatch <- shiny::renderText({
            shiny::req(is(selectedOligoMatch(), "RprimerMatchOligo"))
            c(
                "<b>One mismatch</b><br>",
                selectedOligoMatch()$idOneMismatch[[1]]
            )
        })

        output$twoMismatches <- shiny::renderText({
            shiny::req(is(selectedOligoMatch(), "RprimerMatchOligo"))
            c(
                "<b>Two mismatches</b><br>",
                selectedOligoMatch()$idTwoMismatches[[1]]
            )
        })

        output$threeMismatches <- shiny::renderText({
            shiny::req(is(selectedOligoMatch(), "RprimerMatchOligo"))
            c(
                "<b>Three mismatches</b><br>",
                selectedOligoMatch()$idThreeMismatches[[1]]
            )
        })

        output$fourOrMoreMismatches <- shiny::renderText({
            shiny::req(is(selectedOligoMatch(), "RprimerMatchOligo"))
            c(
                "<b>Four or more mismatches</b><br>",
                selectedOligoMatch()$idFourOrMoreMismatches[[1]]
            )
        })

        list(data = shiny::reactive(ols()))
    })
}
