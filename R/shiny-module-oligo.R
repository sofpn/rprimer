oligoUI <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h5("Design oligos (step 2/5)"),
        shiny::hr(),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                width = 3,
                shiny::h5(shiny::tags$b("Primer settings")),
                shiny::hr(),
                shiny::sliderInput(
                    ns("lengthPrimer"),
                    "Length",
                    value = c(18, 22), min = 15, max = 40
                ),
                shiny::numericInput(
                    ns("maxDegeneracyPrimer"),
                    "Maximum degeneracy (number of sequence variants, 1-64)",
                    value = 4, min = 1, max = 64
                ),
                shiny::checkboxInput(
                    ns("avoidThreeEndRunsPrimer"),
                    "Avoid 3' end runs",
                    value = TRUE
                ),
                shiny::h6(
                    "If primers with more than two runs of the same
                          nucleotide at the terminal 3' end should be avoided
                          (to reduce the risk of mispriming)"
                ),
                shiny::br(),
                shiny::checkboxInput(
                    ns("gcClampPrimer"),
                    "Use GC-clamp",
                    value = TRUE
                ),
                shiny::h6(
                    "A GC-clamp is identified as two to three G or C:s
                    within the last five bases (3' end) of the primer"
                ),
                shiny::br(),
                shiny::sliderInput(
                    ns("gcPrimer"),
                    "GC-content range",
                    value = c(0.4, 0.65), min = 0, max = 1
                ),
                shiny::sliderInput(
                    ns("tmPrimer"),
                    "Melting temperature (Tm) range (Celcius degrees)",
                    value = c(50, 65), min = 20, max = 90
                ),
                shiny::numericInput(
                    ns("concPrimer"),
                    "Primer concentration (20-2000 nM) (for Tm calculation)",
                    value = 500, min = 20, max = 2000
                ),
                shiny::radioButtons(
                    ns("designStrategyPrimer"),
                    "Primer design strategy",
                    choices = c(
                        "Ambiguous" = "ambiguous",
                        "Mixed" = "mixed"
                    ),
                    selected = "ambiguous"
                ),
                shiny::h6(
                    "Ambigous primers can have degenerate bases at any position"
                ),
                shiny::h6(
                    "Mixed primers consist of a shorter degenerate part at the 3' end
                    (1/3 of the primer) and a longer consensus part at the 5'
                    end (2/3 of the primer), which instead of having ambiguous
                    bases contains the most frequently occuring nucleotides"
                ),
                shiny::br(),
                shiny::h5(shiny::tags$b("Probe settings")),
                shiny::hr(),
                shiny::checkboxInput(
                    ns("probe"),
                    "Design probes",
                    value = FALSE
                ),
                shiny::conditionalPanel(
                    ns = shiny::NS(id),
                    condition = "input.probe == true",
                    shiny::sliderInput(
                        ns("lengthProbe"),
                        "Length",
                        value = c(18, 22), min = 15, max = 40
                    ),
                    shiny::numericInput(
                        ns("maxDegeneracyProbe"),
                        "Maximum degeneracy (number of sequence variants, 1-64)",
                        value = 4, min = 1, max = 64
                    ),
                    shiny::checkboxInput(
                        ns("avoidFiveEndGProbe"),
                        "Avoid probes with a 5' end G",
                        value = TRUE
                    ),
                    shiny::sliderInput(
                        ns("gcProbe"),
                        "GC-content range",
                        value = c(0.4, 0.65), min = 0, max = 1,
                    ),
                    shiny::sliderInput(
                        ns("tmProbe"),
                        "Melting temperature range (Celcius degrees)",
                        value = c(50, 70), min = 20, max = 90
                    ),
                    shiny::numericInput(
                        ns("concProbe"),
                        "Probe concentration (20-2000 nM) (for Tm calculation)",
                        value = 250, min = 20, max = 2000
                    )
                ),
                shiny::br(),
                shiny::h5(shiny::tags$b("General settings")),
                shiny::hr(),
                shiny::numericInput(
                    ns("maxGapFrequency"),
                    "Maximum allowed gap proportion at binding sites in
                        target alignment (0-1)",
                    value = 0.01, min = 0, max = 0.2
                ),
                shiny::numericInput(
                    ns("concNa"),
                    "Sodium ion concentration (0.01-1 M) (for Tm calculation)",
                    value = 0.05, min = 0, max = 1
                ),
                shiny::hr(),
                shiny::actionButton(
                    ns("getOligos"), "Get oligos",
                    class = "btn btn-primary"
                )
            ),
            shiny::mainPanel(width = 9,
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
                                     shiny::h6("Select a row for more details"),
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
                                             shiny::h6(shiny::tags$b("Overview")),
                                             DT::dataTableOutput(ns("overviewTable")),
                                             shiny::br(),
                                             shiny::h6(shiny::tags$b("All sequence variants")),
                                             DT::dataTableOutput(ns("allVariantTable")),
                                             shiny::br(),
                                             shiny::h6(shiny::tags$b("Nucleotide distribution in target alignment")),
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
                                             shiny::h6(shiny::tags$b("Proportion of matching sequences")),
                                             shiny::br(),
                                             DT::dataTableOutput(ns("matchTable")),
                                             shiny::br(),
                                             shiny::column(
                                                 width = 12, align = "center",
                                                 spinnerPlot(ns("matchPlot"), width = "75%")
                                             ),
                                             shiny::br(),
                                             shiny::h6(shiny::tags$b("Sequence names")),
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

oligoServer <- function(id, alignment, consensus) {
    shiny::moduleServer(id, function(input, output, session) {
        oligo <- shiny::eventReactive(input$getOligos, {
            shiny::req(consensus())

            validDegenPrimer <- input$maxDegeneracyPrimer >= 1 &&
                input$maxDegeneracyPrimer <= 64
            shinyFeedback::feedbackDanger(
                "maxDegeneracyPrimer", !validDegenPrimer,
                "Enter a value from 1 to 64"
            )

            validDegenProbe <- input$maxDegeneracyProbe >= 1 &&
                input$maxDegeneracyProbe <= 64
            shinyFeedback::feedbackDanger(
                "maxDegeneracyProbe", !validDegenProbe,
                "Enter a value from 1 to 64"
            )

            validConcPrimer <- input$concPrimer >= 20 &&
                input$concPrimer <= 2000
            shinyFeedback::feedbackDanger(
                "concPrimer", !validConcPrimer,
                "Enter a value from 20 to 2000"
            )

            validConcProbe <- input$concProbe >= 20 &&
                input$concProbe <= 2000
            shinyFeedback::feedbackDanger(
                "concProbe", !validConcProbe,
                "Enter a value from 20 to 2000"
            )

            validConcNa <- input$concNa >= 0.01 && input$concNa <= 1
            shinyFeedback::feedbackDanger(
                "concNa", !validConcNa,
                "Enter a value from 0.01 to 1"
            )

            validGap <- input$maxGapFrequency >= 0 &&
                input$maxGapFrequency <= 1
            shinyFeedback::feedbackDanger(
                "maxGapFrequency", !validGap,
                "Enter a value from 0 to 1"
            )

            tryCatch(
                {
                    designOligos(consensus(),
                                 maxGapFrequency = input$maxGapFrequency,
                                 lengthPrimer = input$lengthPrimer,
                                 maxDegeneracyPrimer = input$maxDegeneracyPrimer,
                                 avoidThreeEndRunsPrimer = input$avoidThreeEndRunsPrimer,
                                 gcClampPrimer = input$gcClampPrimer,
                                 gcPrimer = input$gcPrimer,
                                 tmPrimer = input$tmPrimer,
                                 concPrimer = input$concPrimer,
                                 designStrategyPrimer = input$designStrategyPrimer,
                                 probe = input$probe,
                                 lengthProbe = input$lengthProbe,
                                 maxDegeneracyProbe = input$maxDegeneracyProbe,
                                 avoidFiveEndGProbe = input$avoidFiveEndGProbe,
                                 gcProbe = input$gcProbe,
                                 tmProbe = input$tmProbe,
                                 concNa = input$concNa
                    )
                },
                error = function(cond) {
                    if (grepl("must", cond)) {
                        shiny::showNotification(
                            "At least one design setting was invalid. \n
Please review.",
                            type = "error",
                            duration = NULL
                        )
                    } else {
                        shiny::showNotification(
                            "No oligos were found.\n
Try to adjust design settings.",
                            type = "error",
                            duration = NULL
                        )
                    }
                },
                silent = TRUE
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
                searching = FALSE, paging = TRUE,
                pageLength = 100,
                scrollX = TRUE, autoWidth = TRUE,
                ordering = TRUE, scrollY = "1000"
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
