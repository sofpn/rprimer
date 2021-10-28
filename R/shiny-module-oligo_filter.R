oligoFilterUI <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h4("Filter oligos (step 3/5)"),
        shiny::hr(),
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
                        spinnerPlot(
                            ns("oligoPlot"),
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
                        DT::dataTableOutput(ns("oligoTable"))
                    ),
                    shiny::tabPanel(
                        title = "Selection",
                        shiny::br(),
                        shiny::tabsetPanel(
                            type = "pills",
                            shiny::tabPanel(
                                title = "Oligo information",
                                shiny::br(),
                                shiny::uiOutput(ns("getDownloadLinkTxtSel")),
                                shiny::uiOutput(ns("getDownloadLinkFastaSel")),
                                shiny::br(),
                                shiny::h5("Overview"),
                                shiny::hr(),
                                DT::dataTableOutput(ns("overviewTable")),
                                shiny::br(),
                                shiny::h5("All sequence variants"),
                                shiny::hr(),
                                DT::dataTableOutput(ns("allVariantTable")),
                                shiny::br(),
                                shiny::h5("Nucleotide distribution in target alignment"),
                                shiny::hr(),
                                shiny::br(),
                                shiny::column(
                                    width = 12, align = "center",
                                    spinnerPlot(
                                        ns("bindingRegionPlot"),
                                        width = "75%"
                                    )
                                )
                            ),
                            shiny::tabPanel(
                                title = "Match details",
                                shiny::br(),
                                shiny::h5("Proportion of matching sequences"),
                                shiny::br(),
                                DT::dataTableOutput(ns("matchTable")),
                                shiny::br(),
                                shiny::column(
                                    width = 12, align = "center",
                                    spinnerPlot(ns("matchPlot"), width = "75%")
                                ),
                                shiny::br(),
                                shiny::h5("Sequence names"),
                                shiny::hr(),
                                shiny::htmlOutput(ns("matchId"))
                            )
                        )
                    )
                )
            )
        ),
        shiny::hr()
    )
}

oligoFilterServer <- function(id, alignment, consensus, allOligos) {
    shiny::moduleServer(id, function(input, output, session) {
        output$fwd <- shiny::renderUI({
            shiny::req(any(allOligos()$type == "primer" & allOligos()$fwd))

            ns <- session$ns

            list(
                shiny::h5("Forward"),
                shiny::hr(),
                numericInputFrom(allOligos(), ns("fwdRegionFrom")),
                numericInputTo(allOligos(), ns("fwdRegionTo")),
                conservationInput(
                    allOligos(), ns("minOligoIdentityFwd"),
                    direction = "fwd", variable = "identity"
                ),
                conservationInput(
                    allOligos(), ns("minOligoCoverageFwd"),
                    direction = "fwd", variable = "coverage"
                )
            )
        })

        output$rev <- shiny::renderUI({
            shiny::req(any(allOligos()$type == "primer" & allOligos()$rev))

            ns <- session$ns

            list(
                shiny::h5("Reverse"),
                shiny::hr(),
                numericInputFrom(allOligos(), ns("revRegionFrom")),
                numericInputTo(allOligos(), ns("revRegionTo")),
                conservationInput(
                    allOligos(), ns("minOligoIdentityRev"),
                    direction = "rev", variable = "identity"
                ),
                conservationInput(
                    allOligos(), ns("minOligoCoverageRev"),
                    direction = "rev", variable = "coverage"
                )
            )
        })

        output$pr <- shiny::renderUI({
            shiny::req(any(allOligos()$type == "probe"))

            ns <- session$ns

            list(
                shiny::h5("Probe"),
                shiny::hr(),
                numericInputFrom(allOligos(), ns("prRegionFrom")),
                numericInputTo(allOligos(), ns("prRegionTo")),
                conservationInput(
                    allOligos(), ns("minOligoIdentityPr"),
                    type = "probe", variable = "identity"
                ),
                conservationInput(
                    allOligos(), ns("minOligoCoveragePr"),
                    type = "probe", variable = "coverage"
                )
            )
        })

        oligo <- shiny::reactive({
            shiny::req(allOligos())
            filterOligos(allOligos(),
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

        output$oligoPlot <- shiny::renderPlot({
            shiny::req(is(oligo(), "RprimerOligo"))
            plotData(oligo())
        })

        output$getDownloadLinkTxt <- shiny::renderUI({
            shiny::req(is(oligo(), "RprimerOligo"))
            ns <- session$ns
            list(
                shiny::downloadLink(
                    ns("downloadTxt"), "Download table as .txt"
                ),
                shiny::br()
            )
        })

        output$getDownloadLinkFasta <- shiny::renderUI({
            shiny::req(is(oligo(), "RprimerOligo"))
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
                paste0("oligo", "-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                utils::write.table(
                    as.data.frame(oligo()), file,
                    quote = FALSE, sep = "\t",
                    row.names = FALSE
                )
            }
        )

        output$downloadFasta <- shiny::downloadHandler(
            filename <- function() {
                paste0("oligo", "-fasta-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                data <- as(oligo(), "DNAStringSet")
                Biostrings::writeXStringSet(data, file)
            }
        )

        output$oligoTable <- DT::renderDataTable(
            {
                shiny::req(is(oligo(), "RprimerOligo"))
                x <- roundDbls(removeListColumns(as.data.frame(oligo())))
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
                ordering = TRUE, scrollY = "300"
            ),
            rownames = FALSE,
            selection = list(mode = "single")
        )

        shiny::observeEvent(input$oligoTable_rows_selected, {
            Sys.sleep(1)
            shiny::updateTabsetPanel(session, "wizard", selected = "Selection")
        })

        selectedOligo <- shiny::reactive({
            shiny::req(is(oligo(), "RprimerOligo"))
            if (!is.null(input$oligoTable_rows_selected)) {
                oligo()[input$oligoTable_rows_selected, ]
            } else {
                NULL
            }
        })


        match <- shiny::reactive({
            shiny::req(selectedOligo())
            shiny::req(alignment())
            checkMatch(selectedOligo(), alignment())
        })

        output$getDownloadLinkTxtSel <- shiny::renderUI({
            shiny::req(is(selectedOligo(), "RprimerOligo"))
            ns <- session$ns
            list(
                shiny::downloadLink(
                    ns("downloadTxtSel"), "Download table as .txt"
                ),
                shiny::br()
            )
        })

        output$getDownloadLinkFastaSel <- shiny::renderUI({
            shiny::req(is(selectedOligo(), "RprimerOligo"))
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
                paste0("oligo-selection", "-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                utils::write.table(
                    as.data.frame(selectedOligo()), file,
                    quote = FALSE, sep = "\t",
                    row.names = FALSE
                )
            }
        )

        output$downloadFastaSel <- shiny::downloadHandler(
            filename <- function() {
                paste0("oligo-selection", "-fasta-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                data <- as(selectedOligo(), "DNAStringSet")
                Biostrings::writeXStringSet(data, file)
            }
        )

        output$overviewTable <- DT::renderDataTable(
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

        output$allVariantTable <- DT::renderDataTable(
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

        output$bindingRegionPlot <- shiny::renderPlot({
            shiny::req(is(selectedOligo(), "RprimerOligo"))
            shiny::req(consensus())
            from <- selectedOligo()$start
            to <- selectedOligo()$end
            bindingRegion <- consensus()[
                consensus()$position >= from & consensus()$position <= to,
            ]
            plotData(bindingRegion, type = "nucleotide")
        })

        output$matchTable <- DT::renderDataTable(
            {
                shiny::req(match())
                x <- roundDbls(removeListColumns(as.data.frame(match())))
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

        output$matchPlot <- shiny::renderPlot({
            shiny::req(match())
            plotData(match())
        })

        output$matchId <- shiny::renderText({
            shiny::req(match())
            x <- match()
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

        list(data = shiny::reactive(oligo()))
    })
}
