assayUI <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h4("Design assays"),
        shiny::hr(),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::h4("Assay settings"),
                shiny::hr(),
                shiny::sliderInput(
                    ns("length"),
                    shiny::h5("Length"),
                    value = c(60, 120), min = 40, max = 5000,
                ),
                shiny::numericInput(
                    ns("tmDifferencePrimers"),
                    shiny::h5(
                        "Maximum melting temperature
                        difference between primers (Celcius degrees)"
                    ),
                    value = 10, min = 0, max = Inf,
                ),
                shiny::hr(),
                shiny::actionButton(ns("getAssays"), "Get assays")
            ),
            shiny::mainPanel(
                shiny::tabsetPanel(
                    id = ns("wizard"),
                    shiny::tabPanel(
                        title = "Plot",
                        shiny::br(),
                        spinnerPlot(
                            ns("assayPlot"),
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
                        DT::dataTableOutput(ns("assayTable"))
                    ),
                    shiny::tabPanel(
                        title = "Selection",
                        assaySelectionUI(ns("assaySelection"))
                    )
                )
            )
        ),
        shiny::hr()
    )
}

assayServer <- function(id, alignment, consensus, oligo) {
    shiny::moduleServer(id, function(input, output, session) {
        assay <- shiny::eventReactive(input$getAssays, {
            shiny::req(oligo())

            tryCatch(
                {
                    if (is.na(oligo()$length[[1]])) {
                        stop("No assays were found.", .call = FALSE)
                    }

                    assays(oligo(),
                        length = input$length,
                        tmDifferencePrimers = as.numeric(input$tmDifferencePrimers)
                    )
                },
                error = function(cond) {
                    shiny::showNotification(
                        "No assays were found.\n
                        Try to adjust oligo or assay design settings.",
                        type = "error",
                        duration = NULL
                    )
                },
                silent = TRUE
            )
        })

        output$assayPlot <- shiny::renderPlot({
            shiny::req(assay())
            plotData(assay())
        })


        output$getDownloadLinkTxt <- displayDownloadHandlerTxt(
            assay(), session
        )

        output$getDownloadLinkFasta <- displayDownloadHandlerFasta(
            assay(), session
        )

        output$downloadTxt <- downloadHandlerTxt(assay(), "assay")

        output$downloadFasta <- downloadHandlerFasta(assay(), "assay")


        shiny::observeEvent(input$assayTable_rows_selected, {
            Sys.sleep(1)
            shiny::updateTabsetPanel(session, "wizard", selected = "Selection")
        })

        selectedAssay <- shiny::reactive({
            shiny::req(assay())
            if (!is.null(input$assayTable_rows_selected)) {
                assay()[input$assayTable_rows_selected, ]
            } else {
                NULL
            }
        })

        output$assayTable <- assayDataTable(
            assay(),
            selection = list(mode = "single"),
            scrollY = "300", ordering = TRUE
        )

        assaySelectionServer(
            "assaySelection", alignment, consensus, selectedAssay
        )

        list(data = shiny::reactive(assay()))
    })
}

## Module app for testing ======================================================

# assayApp <- function() {
#    data("exampleRprimerAlignment")
#    x <- reactive(exampleRprimerAlignment)
#    data("exampleRprimerProfile")
#    y <- reactive(exampleRprimerProfile)
#    data("exampleRprimerOligo")
#    z <- reactive(exampleRprimerOligo)
#    ui <- fluidPage(
#        assayUI("id")
#    )
#    server <- function(input, output, session) {
#        assayServer("id", alignment = x, consensus = y, oligo = z)
#    }
#    shinyApp(ui, server)
# }
