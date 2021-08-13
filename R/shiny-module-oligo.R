oligoUI <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::titlePanel("Step 2/5: Design oligos"),
        shiny::br(),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::h4("Primer settings"),
                shiny::hr(),
                shiny::sliderInput(
                    ns("lengthPrimer"),
                    shiny::h5("Length"),
                    value = c(18, 22), min = 15, max = 40
                ),
                shiny::numericInput(
                    ns("maxDegeneracyPrimer"),
                    shiny::h5("Maximum degeneracy (1-64)"),
                    value = 4, min = 1, max = 64
                ),
                shiny::checkboxInput(
                    ns("avoidThreeEndRunsPrimer"),
                    shiny::h5("Avoid 3' end runs"),
                    value = TRUE
                ),
                shiny::checkboxInput(
                    ns("gcClampPrimer"),
                    shiny::h5("Use GC clamp"),
                    value = TRUE
                ),
                shiny::sliderInput(
                    ns("gcPrimer"),
                    shiny::h5("GC content range"),
                    value = c(0.4, 0.65), min = 0, max = 1
                ),
                shiny::sliderInput(
                    ns("tmPrimer"),
                    shiny::h5("Melting temperature range (Celcius degrees)"),
                    value = c(50, 65), min = 20, max = 90
                ),
                shiny::numericInput(
                    ns("concPrimer"),
                    shiny::h5("Concentration (20-2000 nM)"),
                    value = 500, min = 20, max = 2000
                ),
                shiny::radioButtons(
                    ns("designStrategyPrimer"),
                    shiny::h5("Design strategy"),
                    choices = c(
                        "Ambiguous" = "ambiguous", "Mixed" = "mixed"
                    ),
                    selected = "ambiguous"
                ),
                shiny::br(),
                shiny::h4("Probe settings"),
                shiny::hr(),
                shiny::checkboxInput(
                    ns("probe"),
                    shiny::h5("Design probes"),
                    value = FALSE
                ),
                shiny::conditionalPanel(
                    ns = shiny::NS(id),
                    condition = "input.probe == true",
                    shiny::sliderInput(
                        ns("lengthProbe"),
                        shiny::h5("Length"),
                        value = c(18, 22), min = 15, max = 40
                    ),
                    shiny::numericInput(
                        ns("maxDegeneracyProbe"),
                        shiny::h5("Maximum degeneracy (1-64)"),
                        value = 4, min = 1, max = 64
                    ),
                    shiny::checkboxInput(
                        ns("avoidFiveEndGProbe"),
                        shiny::h5("Avoid 5' end G"),
                        value = TRUE
                    ),
                    shiny::sliderInput(
                        ns("gcProbe"),
                        shiny::h5("GC content range"),
                        value = c(0.4, 0.65), min = 0, max = 1,
                    ),
                    shiny::sliderInput(
                        ns("tmProbe"),
                        shiny::h5("Melting temperature range (Celcius degrees)"),
                        value = c(50, 70), min = 20, max = 90
                    ),
                    shiny::numericInput(
                        ns("concProbe"),
                        shiny::h5("Concentration (20-2000 nM)"),
                        value = 250, min = 20, max = 2000
                    )
                ),
                shiny::br(),
                shiny::h4("General settings"),
                shiny::hr(),
                shiny::numericInput(
                    ns("maxGapFrequency"),
                    shiny::h5("Maximum gap frequency (0-1)"),
                    value = 0.01, min = 0, max = 0.2
                ),
                shiny::numericInput(
                    ns("concNa"),
                    shiny::h5("Sodium ion concentration (0.01-1 M)"),
                    value = 0.05, min = 0, max = 1
                ),
                shiny::hr(),
                shiny::actionButton(ns("getOligos"), "Get oligos")
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
                                    "Proportion of matching sequences within the
                                                  intended target binding region in the input alignment"
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

oligoServer <- function(id, alignment, consensus) {
    shiny::moduleServer(id, function(input, output, session) {
        ols <- shiny::eventReactive(input$getOligos, {
            shiny::req(is(consensus(), "RprimerProfile"))

            validDegenPrimer <- input$maxDegeneracyPrimer >= 1 && input$maxDegeneracyPrimer <= 64
            shinyFeedback::feedbackDanger(
                "maxDegeneracyPrimer", !validDegenPrimer, "Enter a value from 1 to 64"
            )

            validDegenProbe <- input$maxDegeneracyProbe >= 1 && input$maxDegeneracyProbe <= 64
            shinyFeedback::feedbackDanger(
                "maxDegeneracyProbe", !validDegenProbe, "Enter a value from 1 to 64"
            )

            validConcPrimer <- input$concPrimer >= 20 && input$concPrimer <= 2000
            shinyFeedback::feedbackDanger(
                "concPrimer", !validConcPrimer, "Enter a value from 20 to 2000"
            )

            validConcProbe <- input$concProbe >= 20 && input$concProbe <= 2000
            shinyFeedback::feedbackDanger(
                "concProbe", !validConcProbe, "Enter a value from 20 to 2000"
            )

            validConcNa <- input$concNa >= 0.01 && input$concNa <= 1
            shinyFeedback::feedbackDanger(
                "concNa", !validConcNa, "Enter a value from 0.01 to 1"
            )

            validGap <- input$maxGapFrequency >= 0 && input$maxGapFrequency <= 1
            shinyFeedback::feedbackDanger(
                "maxGapFrequency", !validGap, "Enter a value from 0 to 1"
            )

            tryCatch(
                {
                    oligos(consensus(),
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
            shiny::req(is(ols(), "RprimerOligo"))
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

        output$downloadFasta <- shiny::downloadHandler(
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
                shiny::req(is(ols(), "RprimerOligo"))
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

        output$downloadSelectionFasta <- shiny::downloadHandler(
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
