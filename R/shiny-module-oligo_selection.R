oligoSelectionUI <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::tabsetPanel(
            shiny::tabPanel(
                title = "Oligo information",
                shiny::br(),
                shiny::uiOutput(ns("getDownloadLinkTxt")),
                shiny::uiOutput(ns("getDownloadLinkFasta")),
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
                shiny::h5(
                    "bla bla at the intended target binding site in the alignment)"
                ),
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

oligoSelectionServer <- function(id, alignment, consensus, oligo) {
    shiny::moduleServer(id, function(input, output, session) {
        match <- shiny::reactive({
            shiny::req(oligo())
            shiny::req(alignment())
            checkMatch(oligo(), alignment())
        })

        output$getDownloadLinkTxt <- displayDownloadHandlerTxt(
            oligo(), session
        )

        output$getDownloadLinkFasta <- displayDownloadHandlerFasta(
            oligo(), session
        )

        output$downloadTxt <- downloadHandlerTxt(oligo(), "oligo_selection")

        output$downloadFasta <- downloadHandlerFasta(oligo(), "oligo_selection")

        output$overviewTable <- oligoDataTable(oligo())

        output$allVariantTable <- oligoAllVariantTable(oligo())

        output$bindingRegionPlot <- nucleotidePlot(oligo(), consensus())

        output$matchTable <- oligoMatchTable(match())

        output$matchPlot <- shiny::renderPlot({
            shiny::req(match())
            plotData(match())
        })

        output$matchId <- printMatchId(match())
    })
}

## Module app for testing ======================================================

# oligoSelectionApp <- function() {
#    data("exampleRprimerAlignment")
#    x <- reactive(exampleRprimerAlignment)
#    data("exampleRprimerProfile")
#    y <- reactive(exampleRprimerProfile)
#    data("exampleRprimerOligo")
#    z <- reactive(exampleRprimerOligo[1, ])
#    ui <- fluidPage(
#        oligoSelectionUI("id")
#    )
#    server <- function(input, output, session) {
#        oligoSelectionServer("id", alignment = x, consensus = y, oligo = z)
#   }
#    shinyApp(ui, server)
# }
