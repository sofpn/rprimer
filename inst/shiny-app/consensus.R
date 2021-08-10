consensusUI <- function(id) {
    ns <- NS(id)

    tagList(
        titlePanel("Step 1/5: Make consensus profile"),
        br(),
        sidebarLayout(
            sidebarPanel(
                h5("Region of interest (position)"),
                uiOutput(ns("roiFrom")),
                uiOutput(ns("roiTo")),
                sliderInput(
                    ns("ambiguityThreshold"),
                    h5("Threshold for an ambiguous base"),
                    value = 0.05, min = 0, max = 0.2,
                ),
                br(),
                htmlOutput(ns("ambiguityExplanation")),
                br(),
                actionButton(ns("getConsensusProfile"), "Get consensus profile"),
            ),
            mainPanel(
                tabsetPanel(
                    tabPanel(
                        title = "Plot",
                        br(),
                        shinycssloaders::withSpinner(
                            plotOutput(
                                ns("consensusPlot"),
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
                        DT::dataTableOutput(ns("consensusTable"))
                    )
                )
            )
        )
    )
}

consensusServer <- function(id, alignment) {
    moduleServer(id, function(input, output, session) {
        output$roiFrom <- renderUI({
            req(is(alignment(), "DNAMultipleAlignment"))
            ns <- session$ns
            numericInput(
                ns("roiFrom"), h5("From"),
                min = 1, max = ncol(alignment()) - 1, value = 1,
            )
        })

        output$roiTo <- renderUI({
            req(is(alignment(), "DNAMultipleAlignment"))
            ns <- session$ns
            numericInput(
                ns("roiTo"), h5("To"),
                min = 2, max = ncol(alignment()), value = ncol(alignment())
            )
        })

        output$ambiguityExplanation <- renderText({
            paste(
                "At each position, all bases with an occurence of more than",
                input$ambiguityThreshold * 100, "% will be included when the
                ambiguous base is determined"
            )
        })

        consensus <- eventReactive(input$getConsensusProfile, {
            req(is(alignment(), "DNAMultipleAlignment"))
            consensus <- consensusProfile(alignment(), input$ambiguityThreshold)
            selection <- consensus[
                consensus$position >= as.numeric(input$roiFrom) &
                    consensus$position <= as.numeric(input$roiTo),
            ]
            selection
        })

        output$consensusPlot <- renderPlot({
            req(is(consensus(), "RprimerProfile"))
            plotData(consensus())
        })

        output$getData <- renderUI({
            req(is(consensus(), "RprimerProfile"))
            ns <- session$ns
            downloadLink(ns("download"), "Download table as .txt")
        })

        output$download <- downloadHandler(
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
                req(is(consensus(), "RprimerProfile"))
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

        list(data = reactive(consensus()))
    })
}
