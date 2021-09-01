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
                        oligoSelectionUI(ns("selection"))
                    )
                )
            )
        ),
        shiny::hr()
    )
}

oligoFilterServer <- function(id, alignment, consensus, oligo) {
    stopifnot(is.reactive(oligo))

    shiny::moduleServer(id, function(input, output, session) {
        output$fwd <- shiny::renderUI({
            shiny::req(any(oligo()$type == "primer" & oligo()$fwd))

            ns <- session$ns

            list(
                shiny::h5("Forward"),
                shiny::hr(),
                numericInputFrom(oligo(), ns("fwdRegionFrom")),
                numericInputTo(oligo(), ns("fwdRegionTo")),
                conservationInput(
                    oligo(), ns("minOligoIdentityFwd"),
                    direction = "fwd", variable = "identity"
                ),
                conservationInput(
                    oligo(), ns("minOligoCoverageFwd"),
                    direction = "fwd", variable = "coverage"
                )
            )
        })

        output$rev <- shiny::renderUI({
            shiny::req(any(oligo()$type == "primer" & oligo()$rev))

            ns <- session$ns

            list(
                shiny::h5("Reverse"),
                shiny::hr(),
                numericInputFrom(oligo(), ns("revRegionFrom")),
                numericInputTo(oligo(), ns("revRegionTo")),
                conservationInput(
                    oligo(), ns("minOligoIdentityRev"),
                    direction = "rev", variable = "identity"
                ),
                conservationInput(
                    oligo(), ns("minOligoCoverageRev"),
                    direction = "rev", variable = "coverage"
                )
            )
        })

        output$pr <- shiny::renderUI({
            shiny::req(any(oligo()$type == "probe"))

            ns <- session$ns

            list(
                shiny::h5("Probe"),
                shiny::hr(),
                numericInputFrom(oligo(), ns("prRegionFrom")),
                numericInputTo(oligo(), ns("prRegionTo")),
                conservationInput(
                    oligo(), ns("minOligoIdentityPr"),
                    type = "probe", variable = "identity"
                ),
                conservationInput(
                    oligo(), ns("minOligoCoveragePr"),
                    type = "probe", variable = "coverage"
                )
            )
        })

        filteredOligo <- shiny::reactive({
            shiny::req(oligo())
            filterOligos(oligo(),
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


        shiny::observeEvent(input$oligoTable_rows_selected, {
            Sys.sleep(1)
            shiny::updateTabsetPanel(session, "wizard", selected = "Selection")
        })

        selectedOligo <- shiny::reactive({
            shiny::req(filteredOligo())
            if (!is.null(input$oligoTable_rows_selected)) {
                filteredOligo()[input$oligoTable_rows_selected, ]
            } else {
                NULL
            }
        })

        output$oligoPlot <- shiny::renderPlot({
            shiny::req(filteredOligo())
            plotData(filteredOligo())
        })

        output$getDownloadLinkTxt <- displayDownloadHandlerTxt(
            filteredOligo(), session
        )

        output$getDownloadLinkFasta <- displayDownloadHandlerFasta(
            filteredOligo(), session
        )

        output$downloadTxt <- downloadHandlerTxt(filteredOligo(), "oligo_filtered")

        output$downloadFasta <- downloadHandlerFasta(filteredOligo(), "oligo_filtered")

        output$oligoTable <- oligoDataTable(
            filteredOligo(),
            selection = list(mode = "single"),
            scrollY = "300", ordering = TRUE
        )

        oligoSelectionServer(
            "selection", alignment, consensus, selectedOligo
        )

        list(data = shiny::reactive(filteredOligo()))
    })
}


## Module app for testing ======================================================

 #oligoFilterApp <- function() {
 #   data("exampleRprimerAlignment")
 #   x <- reactive(exampleRprimerAlignment)
 #   data("exampleRprimerProfile")
 #   y <- reactive(exampleRprimerProfile)
 #   data("exampleRprimerOligo")
 #   z <- reactive(exampleRprimerOligo)
 #   ui <- fluidPage(
 #       oligoFilterUI("id")
 #   )
 #   server <- function(input, output, session) {
 #       oligoFilterServer("id", alignment = x, consensus = y, oligo = z)
 #   }
 #   shinyApp(ui, server)
 #}
