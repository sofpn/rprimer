oligoUI <- function(id) {
    ns <- NS(id)

    tagList(
        titlePanel("Step 2/5: Design oligos"),
        br(),
        sidebarLayout(
            sidebarPanel(
                h4("Primer settings"),
                hr(),
                sliderInput(
                    ns("lengthPrimer"),
                    h5("Length"),
                    value = c(18, 22), min = 15, max = 40
                ),
                numericInput(
                    ns("maxDegeneracyPrimer"),
                    h5("Maximum degeneracy (1-64)"),
                    value = 4, min = 1, max = 64
                ),
                checkboxInput(
                    ns("avoidThreeEndRunsPrimer"),
                    h5("Avoid 3' end runs"),
                    value = TRUE
                ),
                checkboxInput(
                    ns("gcClampPrimer"),
                    h5("Use GC clamp"),
                    value = TRUE
                ),
                sliderInput(
                    ns("gcPrimer"),
                    h5("GC content range"),
                    value = c(0.4, 0.65), min = 0, max = 1
                ),
                sliderInput(
                    ns("tmPrimer"),
                    h5("Melting temperature range (Celcius degrees)"),
                    value = c(50, 65), min = 20, max = 90
                ),
                numericInput(
                    ns("concPrimer"),
                    h5("Concentration (20-2000 nM)"),
                    value = 500, min = 20, max = 2000
                ),
                radioButtons(
                    ns("designStrategyPrimer"),
                    h5("Design strategy"),
                    choices = c(
                        "Ambiguous" = "ambiguous", "Mixed" = "mixed"
                    ),
                    selected = "ambiguous"
                ),
                br(),
                h4("Probe settings"),
                hr(),
                checkboxInput(
                    ns("probe"),
                    h5("Design probes"),
                    value = FALSE
                ),
                conditionalPanel(
                    ns = NS(id),
                    condition = "input.probe == true",
                    sliderInput(
                        ns("lengthProbe"),
                        h5("Length"),
                        value = c(18, 22), min = 15, max = 40
                    ),
                    numericInput(
                        ns("maxDegeneracyProbe"),
                        h5("Maximum degeneracy (1-64)"),
                        value = 4, min = 1, max = 64
                    ),
                    checkboxInput(
                        ns("avoidFiveEndGProbe"),
                        h5("Avoid 5' end G"),
                        value = TRUE
                    ),
                    sliderInput(
                        ns("gcProbe"),
                        h5("GC content range"),
                        value = c(0.4, 0.65), min = 0, max = 1,
                    ),
                    sliderInput(
                        ns("tmProbe"),
                        h5("Melting temperature range (Celcius degrees)"),
                        value = c(50, 70), min = 20, max = 90
                    ),
                    numericInput(
                        ns("concProbe"),
                        h5("Concentration (20-2000 nM)"),
                        value = 250, min = 20, max = 2000
                    )
                ),
                br(),
                h4("General settings"),
                hr(),
                numericInput(
                    ns("maxGapFrequency"),
                    h5("Maximum gap frequency (0-1)"),
                    value = 0.01, min = 0, max = 0.2
                ),
                numericInput(
                    ns("concNa"),
                    h5("Sodium ion concentration (0.01-1 M)"),
                    value = 0.05, min = 0, max = 1
                ),
                hr(),
                actionButton(ns("getOligos"), "Get oligos")
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
                                    "Proportion of matching sequences within the
                                                  intended target binding region in the input alignment"
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

oligoServer <- function(id, alignment, consensus) {
    moduleServer(id, function(input, output, session) {
        ols <- eventReactive(input$getOligos, {
            req(is(consensus(), "RprimerProfile"))

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
                        showNotification(
                            "At least one design setting was invalid. \n
                                  Please review.",
                            type = "error",
                            duration = NULL
                        )
                    } else {
                        showNotification(
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
            req(is(ols(), "RprimerOligo"))
            ns <- session$ns
            list(
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

        output$getSelectionData <- renderUI({
            req(is(selectedOligo(), "RprimerOligo"))
            ns <- session$ns
            list(
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
