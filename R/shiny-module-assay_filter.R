assayFilterUI <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h4("Filter assays"),
        shiny::hr(),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::uiOutput(ns("assayFilter"))
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

assayFilterServer <- function(id, alignment, consensus, assay) {
    shiny::moduleServer(id, function(input, output, session) {
        filteredAssay <- shiny::reactive({
            shiny::req(assay())
            x <- assay()
            x <- as.data.frame(x)
            emptyRow <- makeEmptyRow(x)
            x <- x[
                x$start >= input$assayRegionFrom & x$end <= input$assayRegionTo,
            ]
            x <- x[x$score <= input$maxAssayScore, ]
            if (nrow(x) == 0L) {
                x <- emptyRow
            }
            RprimerAssay(x)
        })

        output$assayFilter <- shiny::renderUI({
            ns <- session$ns

            list(
                numericInputFrom(assay(), ns("assayRegionFrom")),
                numericInputTo(assay(), ns("assayRegionTo")),
                sliderInput(ns("maxAssayScore"),
                    shiny::h5("Maximum score (lower is better)"),
                    min = min(assay()$score, na.rm = TRUE),
                    max = max(assay()$score, na.rm = TRUE),
                    value = max(assay()$score, na.rm = TRUE)
                )
            )
        })

        output$assayPlot <- shiny::renderPlot({
            shiny::req(filteredAssay())
            plotData(filteredAssay())
        })


        output$getDownloadLinkTxt <- displayDownloadHandlerTxt(
            filteredAssay(), session
        )

        output$getDownloadLinkFasta <- displayDownloadHandlerFasta(
            filteredAssay(), session
        )

        output$downloadTxt <- downloadHandlerTxt(
            filteredAssay(), "assay_filtered"
        )

        output$downloadFasta <- downloadHandlerFasta(
            filteredAssay(), "assay_filtered"
        )

        shiny::observeEvent(input$assayTable_rows_selected, {
            Sys.sleep(1)
            shiny::updateTabsetPanel(session, "wizard", selected = "Selection")
        })

        selectedAssay <- shiny::reactive({
            shiny::req(filteredAssay())
            if (!is.null(input$assayTable_rows_selected)) {
                filteredAssay()[input$assayTable_rows_selected, ]
            } else {
                NULL
            }
        })

        output$assayTable <- assayDataTable(
            filteredAssay(),
            selection = list(mode = "single"),
            scrollY = "300", ordering = TRUE
        )

        assaySelectionServer(
            "assaySelection", alignment, consensus, selectedAssay
        )

        list(data = shiny::reactive(filteredAssay()))
    })
}

## Module app for testing ======================================================

# assayFilterApp <- function() {
#    data("exampleRprimerAlignment")
#    x <- reactive(exampleRprimerAlignment)
#    data("exampleRprimerProfile")
#    y <- reactive(exampleRprimerProfile)
#    data("exampleRprimerAssay")
#    z <- reactive(exampleRprimerAssay)
#    ui <- fluidPage(
#        assayFilterUI("id")
#    )
#    server <- function(input, output, session) {
#        assayFilterServer("id", alignment = x, consensus = y, assay = z)
#    }
#    shinyApp(ui, server)
# }
