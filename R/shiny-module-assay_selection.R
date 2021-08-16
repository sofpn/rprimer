assaySelectionUI <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::br(),
        shiny::h5("Assay information"),
        br(),
        shiny::uiOutput(ns("getDownloadLinkTxt")),
        shiny::uiOutput(ns("getDownloadLinkFasta")),
        shiny::br(),
        DT::dataTableOutput(ns("overviewTable")),
        shiny::br(),
        shiny::uiOutput(ns("detailsTab"))
    )
}

assaySelectionServer <- function(id, alignment, consensus, assay) {
    shiny::moduleServer(id, function(input, output, session) {
        assayList <- shiny::reactive({
            shiny::req(assay())
            splitAssayToList(assay())
        })

        assayMatch <- shiny::reactive({
            shiny::req(assay())
            shiny::req(alignment())
            if (is.na(assay()$length[[1]])) {
                NULL
            } else {
                checkMatch(assay(), alignment())
            }
        })

        assayMatchList <- shiny::reactive({
            shiny::req(assayMatch())
            splitAssayToList(assayMatch())
        })

        output$getDownloadLinkTxt <- displayDownloadHandlerTxt(
            assay(), session
        )

        output$getDownloadLinkFasta <- displayDownloadHandlerFasta(
            assay(), session
        )

        output$downloadTxt <- downloadHandlerTxt(assay(), "assay_selection")

        output$downloadFasta <- downloadHandlerFasta(assay(), "assay_selection")

        output$overviewTable <- assayOverviewTable(assay())

        output$detailsTab <- shiny::renderUI({
            ns <- session$ns
            tabs <- list(
                shiny::tabPanel(
                    title = "Forward",
                    br(),
                    oligoSelectionUI(ns("fwd"))
                ),
                shiny::tabPanel(
                    title = "Reverse",
                    br(),
                    oligoSelectionUI(ns("rev"))
                ),
                shiny::tabPanel(
                    title = "Probe",
                    br(),
                    oligoSelectionUI(ns("pr"))
                )
            )

            if (length(assayList()) == 3) {
                do.call(shiny::tabsetPanel, tabs)
            } else {
                do.call(shiny::tabsetPanel, tabs[seq_len(2)])
            }
        })

        oligoFwd <- reactive(convertToOligo(assayList()[[1]]))

        oligoRev <- reactive(convertToOligo(assayList()[[2]]))

        oligoPr <- reactive(
            if (length(assayList()) == 3) {
                convertToOligo(assayList()[[3]])
            } else {
                NULL
            }
        )

        oligoSelectionServer("fwd", alignment, consensus, oligoFwd)

        oligoSelectionServer("rev", alignment, consensus, oligoRev)

        oligoSelectionServer("pr", alignment, consensus, oligoPr)
    })
}

## Module app for testing ======================================================

assaySelectionApp <- function() {
    data("exampleRprimerAlignment")
    x <- reactive(exampleRprimerAlignment)
    data("exampleRprimerProfile")
    y <- reactive(exampleRprimerProfile)
    data("exampleRprimerAssay")
    z <- reactive(exampleRprimerAssay[1, ])
    ui <- fluidPage(
        assaySelectionUI("id")
    )
    server <- function(input, output, session) {
        assaySelectionServer("id", alignment = x, consensus = y, assay = z)
    }
    shinyApp(ui, server)
}
