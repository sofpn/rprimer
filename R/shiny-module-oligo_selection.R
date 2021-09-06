oligoSelectionUI <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::tabsetPanel(
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
}

oligoSelectionServer <- function(id, alignment, consensus, selectedOligo) {
    shiny::moduleServer(id, function(input, output, session) {
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
                if (is.na(x$Score[[1]])) {
                    x <- x[!names(x) %in% "Score"]
                }
                x
            },
            options = list(
                info = FALSE,
                searching = FALSE, paging = FALSE,
                scrollX = TRUE, autoWidth = TRUE,
                ordering = TRUE
            ),
            rownames = FALSE,
            selection = "none"
        )

        output$allVariantTable <- DT::renderDataTable(
            {
                shiny::req(is(selectedOligo(), "RprimerOligo"))
                x <- roundDbls(makeListTable(as.data.frame(selectedOligo())))
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
    })
}
