consensusUI <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h4("Make consensus profile"),
        shiny::hr(),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::h5("Region of interest (position)"),
                shiny::uiOutput(ns("roiFrom")),
                shiny::uiOutput(ns("roiTo")),
                shiny::sliderInput(
                    ns("ambiguityThreshold"),
                    shiny::h5("Threshold for an ambiguous base"),
                    value = 0.05, min = 0, max = 0.2,
                ),
                shiny::br(),
                shiny::htmlOutput(ns("ambiguityExplanation")),
                shiny::br(),
                shiny::actionButton(
                    ns("getConsensusProfile"), "Get consensus profile"
                ),
            ),
            shiny::mainPanel(
                shiny::tabsetPanel(
                    shiny::tabPanel(
                        title = "Plot",
                        shiny::br(),
                        spinnerPlot(ns("consensusPlot"))
                    ),
                    shiny::tabPanel(
                        title = "Table",
                        shiny::br(),
                        shiny::uiOutput(ns("getDownloadLink")),
                        shiny::br(),
                        DT::dataTableOutput(ns("consensusTable"))
                    )
                )
            )
        ),
        shiny::hr()
    )
}

consensusServer <- function(id, alignment) {
    shiny::moduleServer(id, function(input, output, session) {
        output$roiFrom <- shiny::renderUI({
            shiny::req(is(alignment(), "DNAMultipleAlignment"))
            ns <- session$ns
            shiny::numericInput(
                ns("roiFrom"), shiny::h5("From"),
                min = 1, max = ncol(alignment()) - 1, value = 1,
            )
        })

        output$roiTo <- shiny::renderUI({
            req(is(alignment(), "DNAMultipleAlignment"))
            ns <- session$ns
            shiny::numericInput(
                ns("roiTo"), shiny::h5("To"),
                min = 2, max = ncol(alignment()), value = ncol(alignment())
            )
        })

        output$ambiguityExplanation <- shiny::renderText({
            paste(
                "At each position, all bases with an occurence of more than",
                input$ambiguityThreshold * 100, "% will be included when the
                ambiguous base is determined"
            )
        })

        consensus <- shiny::eventReactive(input$getConsensusProfile, {
            shiny::req(is(alignment(), "DNAMultipleAlignment"))
            consensus <- consensusProfile(alignment(), input$ambiguityThreshold)
            consensus[
                consensus$position >= as.numeric(input$roiFrom) &
                    consensus$position <= as.numeric(input$roiTo),
            ]
        })

        output$consensusPlot <- shiny::renderPlot({
            shiny::req(is(consensus(), "RprimerProfile"))
            plotData(consensus())
        })

        output$getDownloadLink <- displayDownloadHandlerTxt(
            consensus(), session
        )

        output$downloadTxt <- downloadHandlerTxt(
            consensus(), "consensus_profile"
        )

        output$consensusTable <- consensusDataTable(consensus())

        list(data = shiny::reactive(consensus()))
    })
}

## Module app for testing ======================================================

# consensusApp <- function() {
#    data("exampleRprimerAlignment")
#    x <- reactive(exampleRprimerAlignment)
#    ui <- fluidPage(
#        consensusUI("id")
#    )
#    server <- function(input, output, session) {
#        consensusServer("id", alignment = x)
#    }
#    shinyApp(ui, server)
# }
