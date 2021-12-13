consensusUI <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h5("Generate consensus profile (step 1/5)"),
        shiny::hr(),
        shiny::sidebarLayout(
            shiny::sidebarPanel(width = 3,
                shiny::h6(shiny::tags$b("Region of interest")),
                shiny::h6("Specify start and end position of the consensus profile"),
                shiny::uiOutput(ns("roiFrom")),
                shiny::uiOutput(ns("roiTo")),
                shiny::sliderInput(
                    ns("ambiguityThreshold"),
                    shiny::h6(shiny::tags$b("Threshold for ambiguous (degenerate) bases")),
                    value = 0.05, min = 0, max = 0.2,
                ),
                shiny::h6(
                    "The threshold specifies the 'detection level' for
                    ambiguous bases. All DNA bases that occur with a
                    proportion higher than the specified value will be
                    included in the IUPAC consensus character"
                ),
                shiny::h6(shiny::tags$b("Apply mask(s)")),
                shiny::h6(
                    "Masked positions will not appear as oligo binding sites.
                No mask will be applied if 'From' and 'To' are set to 0"
                ),
                shiny::h6("Mask 1"),
                shiny::numericInput(ns("mask1From"), label = "From", value = 0),
                shiny::numericInput(ns("mask1To"), label = "To",  value = 0),
                shiny::h6("Mask 2"),
                shiny::numericInput(ns("mask2From"), label = "From", value = 0),
                shiny::numericInput(ns("mask2To"), label = "To", value = 0),
                shiny::h6("Mask 3"),
                shiny::numericInput(ns("mask3From"), label = "From",  value = 0),
                shiny::numericInput(ns("mask3To"), label = "To", value = 0),
                shiny::br(),
                shiny::actionButton(
                    ns("getConsensusProfile"), "Get consensus profile",
                    class = "btn btn-primary"
                ),
            ),
            shiny::mainPanel(width = 9,
                shiny::tabsetPanel(
                    shiny::tabPanel(
                        title = "Plot",
                        shiny::br(),
                        spinnerPlot(ns("consensusPlot"))
                    ),
                    shiny::tabPanel(
                        title = "Table",
                        shiny::br(),
                        shiny::uiOutput(ns("getDownloadLinkTxt")),
                        shiny::br(),
                        DT::dataTableOutput(ns("consensusTable"))
                    )
                )
            )
        ),
        shiny::hr()
    )
}

consensusServer <- function(id, inputAlignment) {
    shiny::moduleServer(id, function(input, output, session) {
        output$roiFrom <- shiny::renderUI({
            shiny::req(inputAlignment())
            ns <- session$ns
            shiny::numericInput(
                ns("roiFrom"), label = "From",
                min = 1, max = ncol(inputAlignment()) - 1, value = 1,
            )
        })

        output$roiTo <- shiny::renderUI({
            req(inputAlignment())
            ns <- session$ns
            shiny::numericInput(
                ns("roiTo"), label = "To",
                min = 2, max = ncol(inputAlignment()), value = ncol(inputAlignment())
            )
        })

        alignment <- shiny::reactive({
            shiny::req(inputAlignment())
            maskedAln <- inputAlignment()
            mask1 <- arrangeMask(input$mask1From, input$mask1To, inputAlignment())
            mask2 <- arrangeMask(input$mask2From, input$mask2To, inputAlignment())
            mask3 <- arrangeMask(input$mask3From, input$mask3To, inputAlignment())
            mask <- unique(c(mask1, mask2, mask3))
            if (!is.null(mask)) {
                Biostrings::colmask(maskedAln) <- mask
            }
            maskedAln
        })

        consensus <- shiny::eventReactive(input$getConsensusProfile, {
            shiny::req(alignment())
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

        output$getDownloadLinkTxt <- shiny::renderUI({
            shiny::req(is(consensus(), "RprimerProfile"))
            ns <- session$ns
            list(
                shiny::downloadLink(
                    ns("downloadTxt"), "Download table as .txt"
                ),
                shiny::br()
            )
        })

        output$downloadTxt <- shiny::downloadHandler(
            filename <- function() {
                paste0("consensus", "-", Sys.Date(), ".txt")
            },
            content <- function(file) {
                utils::write.table(
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
                    "Identity", "IUPAC", "Coverage"
                )
                x
            },
            options = list(
                info = FALSE,
                searching = FALSE, paging = TRUE,
                pageLength = 100,
                scrollX = TRUE, autoWidth = FALSE,
                ordering = FALSE, scrollY = "1000"
            ),
            rownames = FALSE,
            selection = "none"
        )

        list(data = shiny::reactive(consensus()))
    })
}
