consensusUI <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::titlePanel("Step 1/5: Make consensus profile"),
        shiny::br(),
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
                        shinycssloaders::withSpinner(
                            shiny::plotOutput(
                                ns("consensusPlot"),
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
                        DT::dataTableOutput(ns("consensusTable"))
                    )
                )
            )
        )
    )
}

consensusServer <- function(id, alignment) {
    shiny::moduleServer(id, function(input, output, session) {
        output$roiFrom <-  shiny::renderUI({
            shiny::req(is(alignment(), "DNAMultipleAlignment"))
            ns <- session$ns
            shiny::numericInput(
                ns("roiFrom"),  shiny::h5("From"),
                min = 1, max = ncol(alignment()) - 1, value = 1,
            )
        })

        output$roiTo <-  shiny::renderUI({
            req(is(alignment(), "DNAMultipleAlignment"))
            ns <- session$ns
            shiny::numericInput(
                ns("roiTo"),  shiny::h5("To"),
                min = 2, max = ncol(alignment()), value = ncol(alignment())
            )
        })

        output$ambiguityExplanation <-  shiny::renderText({
            paste(
                "At each position, all bases with an occurence of more than",
                input$ambiguityThreshold * 100, "% will be included when the
                ambiguous base is determined"
            )
        })

        consensus <-  shiny::eventReactive(input$getConsensusProfile, {
            shiny::req(is(alignment(), "DNAMultipleAlignment"))
            consensus <- consensusProfile(alignment(), input$ambiguityThreshold)
            consensus[
                consensus$position >= as.numeric(input$roiFrom) &
                    consensus$position <= as.numeric(input$roiTo),
            ]
        })

        output$consensusPlot <-  shiny::renderPlot({
            shiny::req(is(consensus(), "RprimerProfile"))
            plotData(consensus())
        })

        output$getData <-  shiny::renderUI({
            shiny::req(is(consensus(), "RprimerProfile"))
            ns <- session$ns
            shiny::downloadLink(ns("download"), "Download table as .txt")
        })

        output$download <-  shiny::downloadHandler(
            filename <- function() {
                paste0("consensus_profile-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                write.table(
                    as.data.frame(consensus()), file,
                    quote = FALSE, sep = "\t",
                    row.names = FALSE
                )
            }
        )

        output$consensusTable <- DT::renderDataTable(
            {
                shiny::req(is(consensus(), "RprimerProfile"))
                x <- roundDbls(as.data.frame(consensus()))
                names(x) <- c(
                    "Position", "A", "C", "G", "T", "Other", "Gaps", "Majority",
                    "Identity", "IUPAC", "Entropy", "Coverage"
                )
                x
            },
            options = list(
                info = FALSE,
                searching = FALSE, paging = TRUE,
                scrollX = TRUE, autoWidth = FALSE,
                ordering = FALSE, scrollY = "300"
            ),
            rownames = FALSE,
            selection = "none"
        )

        list(data = shiny::reactive(consensus()))
    })
}
