oligoUI <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h4("Design oligos (step 2/5)"),
        shiny::hr(),
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
                    shiny::h5(
                        "Maximum allowed gap proportion at binding sites in
                        target alignment (0-1)"
                    ),
                    value = 0.01, min = 0, max = 0.2
                ),
                shiny::numericInput(
                    ns("concNa"),
                    shiny::h5("Sodium ion concentration (0.01-1 M)"),
                    value = 0.05, min = 0, max = 1
                ),
                shiny::hr(),
                shiny::actionButton(
                    ns("getOligos"), "Get oligos",
                    class = "btn btn-primary"
                )
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
                        oligoSelectionUI(ns("selection"))
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

        shiny::observeEvent(input$oligoTable_rows_selected, {
            Sys.sleep(1)
            shiny::updateTabsetPanel(session, "wizard", selected = "Selection")
        })

        selectedOligo <- shiny::reactive({
            shiny::req(oligo())
            if (!is.null(input$oligoTable_rows_selected)) {
                oligo()[input$oligoTable_rows_selected, ]
            } else {
                NULL
            }
        })

        output$oligoPlot <- shiny::renderPlot({
            shiny::req(oligo())
            plotData(oligo())
        })

        output$getDownloadLinkTxt <- displayDownloadHandlerTxt(
            oligo(), session
        )

        output$getDownloadLinkFasta <- displayDownloadHandlerFasta(
            oligo(), session
        )

        output$downloadTxt <- downloadHandlerTxt(oligo(), "oligo")

        output$downloadFasta <- downloadHandlerFasta(oligo(), "oligo")

        output$oligoTable <- oligoDataTable(
            oligo(),
            selection = list(mode = "single"),
            scrollY = "300", ordering = TRUE
        )

        oligoSelectionServer(id, alignment, consensus, shiny::reactive(selectedOligo()))

        list(data = shiny::reactive(oligo())
        )
    })
}


## Module app for testing ======================================================

# oligoApp <- function() {
#    data("exampleRprimerAlignment")
#    x <- reactive(exampleRprimerAlignment)
#    data("exampleRprimerProfile")
#    y <- reactive(exampleRprimerProfile)
#    ui <- fluidPage(
#        oligoUI("id")
#    )
#    server <- function(input, output, session) {
#
#        oligo <- oligoServer("id", alignment = x, consensus = y)
#        oligoSelectionServer(
#            "selection", alignment = x, consensus = y, oligo = oligo$selectedOligo
#        )
#
#    }
#    shinyApp(ui, server)
# }
