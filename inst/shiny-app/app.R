library(shiny)
library(shinydashboard)
library(shinycssloaders)
devtools::load_all("C:/Users/Sofia/Desktop/rprimer")

data("exampleRprimerAlignment")

plotWidth = 800
plotHeight = 600

# Custom functions =============================================================

roundDbls <- function(x) {
    x <- as.data.frame(x)
    y <- lapply(seq_len(ncol(x)), function(i) {
        out <- if (is.double(x[, i])) {
            round(x[, i], 2)
        } else if (is.list(x[, i])) {
            if (is.double(x[, i][[1]])) {
                list(round(x[, i][[1]], 2))
            } else {
                list(x[, i][[1]])
            }
        } else {
            x[, i]
        }
    })
    names(y) <- names(x)
    data.frame(do.call("cbind", y))
}

removeListColumns <- function(x) {
    lists <- vapply(seq_len(ncol(x)), function(i) {
        is.list(x[, i])
    }, logical(1L))
    x[, !lists]
}

# Conditional panels ===========================================================

upload <- conditionalPanel(
    condition = "input.fileInput == 'Upload file'",
    br(),
    h5("Upload a multiple DNA sequence alignment in fasta format"),
    fileInput(
        "file1", "", multiple = FALSE, accept = c("text")
    )
)

useExample <- conditionalPanel(
    condition = "input.fileInput == 'Use example data'",
    br(),
    h5("Hepatitis E virus example data")
)

useProbe <- conditionalPanel(
    condition = "input.probe == true",
    sliderInput(
        "lengthProbe",
        h5("Length"),
        value = c(18, 22), min = 15, max = 40,
        width = 200
    ),
    numericInput(
        "maxDegeneracyProbe",
        h5("Maximum degeneracy"),
        value = 4, min = 1, max = 64,
        width = 200
    ),
    checkboxInput(
        "avoidFiveEndGProbe",
        h5("Avoid 5' end G"),
        value = TRUE,
        width = 200
    ),
    sliderInput(
        "gcRangeProbe",
        h5("GC content range"),
        value = c(0.4, 0.65), min = 0, max = 1,
        width = 200
    ),
    sliderInput(
        "tmRangeProbe",
        h5("Melting temperature range"),
        value = c(50, 70), min = 20, max = 90,
        width = 200
    ),
    numericInput(
        "concProbe",
        h5("Concentration (nM)"),
        value = 250, min = 20, max = 2000,
        width = 200
    )
)

# UI ===========================================================================

ui <- dashboardPage(
    dashboardHeader(title = ""),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Target alignment", tabName = "import"),
            menuItem("1.  Consensus profile", tabName = "consensus"),
            menuItem("2.  Oligos", tabName = "oligos"),
            menuItem("3.  Assays", tabName = "assays")
        )
    ),
    dashboardBody(
        tabItems(

        #### Target alignment ####

        tabItem(tabName = "import",
            fluidRow(
            box(width = 12, title = "File import",
                hr(),
                radioButtons(
                    "fileInput", label = NULL,
                    choices = c(
                        "Upload file", "Use example data"
                    ),
                    selected = "Upload file"
                ),
                hr(),
                upload,
                useExample,
                hr(), ## UI output sequence names
                htmlOutput("html1"),
                htmlOutput("html2")
            )
            )
        ),

        #### Consensus #####

        tabItem(tabName = "consensus",
            fluidRow(
            box(width = 12, title = "Consensus profile",
                hr(),
                column(width = 2,
                h4("Settings"),
                br(),
                numericInput(
                    "ambiguityThreshold",
                    h5("Threshold for an ambiguous base"),
                    value = 0.05, min = 0, max = 0.2, width = 200
                ),
                br(),
                uiOutput("roi"),
                uiOutput("roiFrom"),
                uiOutput("roiTo"),
                br(),
                uiOutput("getConsensusProfile")
                ),
                column(width = 10,
                h4("Output"),
                br(),
                uiOutput("showRegion"),
                column(width = 6,
                       column(width = 6,
                uiOutput("zoomFrom")),
                column(width = 6,
                uiOutput("zoomTo"))
                ),
                column(width = 6),
                br(),
                br(),
                tabBox(title = "",
                       width = 12,
                       tabPanel(title = "Plot",
                                br(),
                                withSpinner(plotOutput(
                                    "plot1",
                                    height = plotHeight,
                                    width = "100%"
                                ), color = "grey"),
                                br(),
                                br(),
                                br(),
                                withSpinner(plotOutput(
                                    "plot2",
                                    height = plotHeight * 0.7,
                                    width = "80%"
                                ), color = NA)
                       ),
                       tabPanel(title = "Table",
                                br(),
                                downloadLink("download1", "Download"),
                                br(),
                                br(),
                                DT::dataTableOutput("table1")
                       )
                )
                )
                )
            )
        ),

        #### Oligos ####

        tabItem(tabName = "oligos",
                fluidRow(
                    box(width = 12, title = "Oligos",
                        hr(),
                        column(width = 2,
                               h4("Settings"),
                               br(),
                               h4("Primers"),
                               hr(),
                               sliderInput(
                                   "lengthPrimer",
                                   h5("Length"),
                                   value = c(18, 22), min = 15, max = 40,
                                   width = 200
                               ),
                               numericInput(
                                   "maxDegeneracyPrimer",
                                   h5("Maximum degeneracy"),
                                   value = 4, min = 1, max = 64,
                                   width = 200
                               ),
                               checkboxInput(
                                   "avoidThreeEndRunsPrimer",
                                   h5("Avoid 3' end runs"),
                                   value = TRUE,
                                   width = 200
                               ),
                               checkboxInput(
                                   "gcClampPrimer",
                                   h5("Use GC clamp"),
                                   value = TRUE,
                                   width = 200
                               ),
                               sliderInput(
                                   "gcRangePrimer",
                                   h5("GC content range"),
                                   value = c(0.4, 0.65), min = 0, max = 1,
                                   width = 200
                               ),
                               sliderInput(
                                   "tmRangePrimer",
                                   h5("Melting temperature range"),
                                   value = c(50, 65), min = 20, max = 90,
                                   width = 200
                               ),
                               numericInput(
                                   "concPrimer",
                                   h5("Concentration (nM)"),
                                   value = 500, min = 20, max = 2000,
                                   width = 200
                               ),
                               radioButtons(
                                   "designStrategyPrimer",
                                   h5("Design strategy"),
                                   choices = c("ambiguous", "mixed"),
                                   selected = "ambiguous"
                               ),
                               br(),
                               h4("Probes"),
                               hr(),
                               checkboxInput(
                                   "probe",
                                   h5("Design probes"),
                                   value = FALSE,
                                   width = 200
                               ),
                               useProbe,
                               br(),
                               h4("General"),
                               hr(),
                               numericInput(
                                   "maxGapFrequency",
                                   h5("Maximum gap frequency"),
                                   value = 0.01, min = 0, max = 0.2,
                                   width = 200
                               ),
                               numericInput(
                                   "concNa",
                                   h5("Sodium ion concentration (M)"),
                                   value = 0.05, min = 0, max = 1,
                                   width = 200
                               ),
                               hr(),
                               uiOutput("getOligos")
                        ),
                        column(width = 10,
                               h4("Output"),
                               br(),
                               uiOutput("showRegionOligos"),
                               column(width = 6,
                               column(width = 6,
                                      uiOutput("fwdRegionFrom"),
                                      uiOutput("revRegionFrom"),
                                      uiOutput("prRegionFrom")

                               ),
                               column(width = 6,
                                      uiOutput("fwdRegionTo"),
                                      uiOutput("revRegionTo"),
                                      uiOutput("prRegionTo"),
                               )
                               ),
                               br(),
                               tabBox(title = "", width = 12,
                                      tabPanel(title = "Plot",
                                               br(),
                                               withSpinner(plotOutput(
                                                   "plot3",
                                                   height = plotHeight,
                                                   width = "100%"
                                               ), color = "grey")
                                      ),
                                      tabPanel(title = "Table and selection",
                                               br(),
                                               uiOutput("allOligos"),
                                               br(),
                                               downloadLink(
                                                   "download2", "Download"
                                               ),
                                               br(),
                                               br(),
                                               DT::dataTableOutput("table2"),
                                               hr(),
                                               br(),
                                               uiOutput("selectedOligo"),
                                               br(),
                                               br(),
                                               withSpinner(plotOutput(
                                                   "plot4",
                                                   height = plotHeight / 2,
                                                   width = "100%"
                                               ), color = "grey"),
                                               br(),
                                               withSpinner(plotOutput(
                                                   "plot5",
                                                   height = plotHeight / 4,
                                                   width = "100%"
                                               ), color = "grey"),
                                               br(),
                                               DT::dataTableOutput("table3"),
                                               br(),
                                               DT::dataTableOutput("table4"),
                                               br(),
                                               htmlOutput("html3")
                                      )
                               )
                        )
                    )
                )
        ),

        #### Assays ####

        tabItem(tabName = "assays",
            box(width = 12, title = "Assays",
                h5("Design settings"),
                actionButton("getAssays", "Design assays"),
                hr()
                #plotOutput("plot4", height = plotHeight, width = plotWidth)
            )
        )



            )))

# Server =======================================================================

server <- function(input, output) {

    #### Data ####

    aln <- reactive({
        if (input$fileInput == "Upload file") {
            req(input$file1)
            x <- Biostrings::readDNAMultipleAlignment(
                input$file1$datapath, format = "fasta"
            )
        } else if (input$fileInput == "Use example data") {
            data("exampleRprimerAlignment")
            x <- exampleRprimerAlignment
        }
        x
    })

    consensus <- eventReactive(input$getConsensusProfile, {
        consensus <- consensusProfile(aln(), input$ambiguityThreshold)
        consensus[
            consensus$position >= as.numeric(input$roiFrom) &
                consensus$position <= as.numeric(input$roiTo),
        ]
    })

    consensusSelection <- reactive({
        consensus()[
            consensus()$position >= as.numeric(input$zoomFrom) &
                consensus()$position <= as.numeric(input$zoomTo),
        ]
    })

    oligoCandidates <- eventReactive(input$getOligos, {
        req(consensus())
        oligos(consensus(),
               maxGapFrequency = input$maxGapFrequency,
               lengthPrimer = input$lengthPrimer,
               maxDegeneracyPrimer = input$maxDegeneracyPrimer,
               avoidThreeEndRunsPrimer = input$avoidThreeEndRunsPrimer,
               gcRangePrimer = input$gcRangePrimer,
               tmRangePrimer = input$tmRangePrimer,
               concPrimer = input$concPrimer,
               designStrategyPrimer = input$designStrategyPrimer,
               probe = input$probe,
               lengthProbe = input$lengthProbe,
               maxDegeneracyProbe = input$maxDegeneracyProbe,
               avoidFiveEndGProbe = input$avoidFiveEndGProbe,
               gcRangeProbe = input$gcRangeProbe,
               tmRangeProbe = input$tmRangeProbe,
               concNa = input$concNa
               )
    })

    oligoSelection <- reactive({
        req(oligoCandidates())
        x <- oligoCandidates()
        fwd <- x$fwd & x$type == "primer" &
            x$start >= as.numeric(input$fwdRegionFrom)  &
            x$end <= as.numeric(input$fwdRegionTo)
        rev <- (x$rev & x$type == "primer" &
            x$start >= as.numeric(input$revRegionFrom)  &
            x$end <= as.numeric(input$revRegionTo))

        x[fwd | rev, ]
      #### If nrow 0 ##################### # fix selection , add probe ###################
    })

    selectedOligo <- reactive({
        req(oligoCandidates())
        oligoCandidates()[input$table2_rows_selected, ]
    })

    selectedOligoMatch <- reactive({
        req(selectedOligo())
        req(aln())
        checkMatch(selectedOligo(), aln())
    })

    assayCandidates <- eventReactive(input$getAssays, {
        req(oligoCandidates())
        assays(oligoCandidates())
    })

    #### Render UI ####

    output$getConsensusProfile <- renderUI({
        req(aln())
        actionButton("getConsensusProfile", "Get consensus profile")
    })

    output$roi <- renderUI({
        req(aln())
        h5("Region of interest")
    })

    output$roiFrom <- renderUI({
        req(aln())
        numericInput(
            "roiFrom", h5("From"),
            min = 1, max = ncol(aln()) - 1, value = 1,
            width = 200
        )
    })

    output$roiTo <- renderUI({
        req(aln())
        numericInput(
            "roiTo", h5("To"), min = 2, max = ncol(aln()), value = ncol(aln()),
            width = 200
        )
    })

    output$showRegion <- renderUI({
        req(consensus())
        h5("Show region")
    })

    output$zoomFrom <- renderUI({
        req(consensus())
        numericInput(
            "zoomFrom", h5("From"),
            min = min(consensus()$position, na.rm = TRUE),
            max = max(consensus()$position, na.rm = TRUE),
            value = min(consensus()$position, na.rm = TRUE)
        )
    })

    output$zoomTo <- renderUI({
        req(consensus())
        numericInput(
            "zoomTo", h5("To"),
            min = min(consensus()$position, na.rm = TRUE),
            max = max(consensus()$position, na.rm = TRUE),
            value = max(consensus()$position, na.rm = TRUE)
        )
    })

    output$getOligos <- renderUI({
        req(consensus())
        actionButton("getOligos", "Get oligos")
    })

    output$showRegionOligos <- renderUI({
        req(oligoCandidates())
        h5("Show region")
    })

    output$fwdRegionFrom <- renderUI({
        req(oligoCandidates())
        numericInput("fwdRegionFrom", h5("Fwd, from"),
            min = oligoCandidates()$roiStart[[1]],
            max = oligoCandidates()$roiEnd[[1]],
            value = oligoCandidates()$roiStart[[1]]
        )
    })

    output$fwdRegionTo <- renderUI({
        req(oligoCandidates())
        numericInput("fwdRegionTo", h5("Fwd, to"),
                     min = oligoCandidates()$roiStart[[1]],
                     max = oligoCandidates()$roiEnd[[1]],
                     value = oligoCandidates()$roiEnd[[1]]
        )
    })


    output$revRegionFrom <- renderUI({
        req(oligoCandidates())
        numericInput("revRegionFrom", h5("Rev, from"),
                     min = oligoCandidates()$roiStart[[1]],
                     max = oligoCandidates()$roiEnd[[1]],
                     value = oligoCandidates()$roiStart[[1]]
        )
    })

    output$revRegionTo <- renderUI({
        req(oligoCandidates())
        numericInput("revRegionTo", h5("Rev, to"),
                     min = oligoCandidates()$roiStart[[1]],
                     max = oligoCandidates()$roiEnd[[1]],
                     value = oligoCandidates()$roiEnd[[1]]
        )
    })

    output$prRegionFrom <- renderUI({
        req(oligoCandidates())
        numericInput("prRegionFrom", h5("Probe, from"),
                     min = oligoCandidates()$roiStart[[1]],
                     max = oligoCandidates()$roiEnd[[1]],
                     value = oligoCandidates()$roiStart[[1]]
        )
    })

    output$prRegionTo <- renderUI({
        req(oligoCandidates())
        numericInput("prRegionTo", h5("Probe, to"),
                     min = oligoCandidates()$roiStart[[1]],
                     max = oligoCandidates()$roiEnd[[1]],
                     value = oligoCandidates()$roiEnd[[1]]
        )
    })

    output$allOligos <- renderUI({
        req(!is.null(oligoCandidates()))
        h4("All oligos")
    })

    output$selectedOligo <- renderUI({
        req(selectedOligo())
        h4("Selection")
    })

    #### Html ####

    output$html1 <- renderText({
        req(aln())
        paste("Number of sequences:", nrow(aln()))
    })

    output$html2 <- renderText({
        req(aln())
        paste("Alignment length:", ncol(aln()))
    })

    output$html3 <- renderText(({
        req(selectedOligoMatch())
        paste("Perfect match:", selectedOligoMatch()$idPerfectMatch)
    }))

    #### Plots ####

    output$plot1 <- renderPlot({
        req(consensusSelection())
        plotData(consensusSelection())
    })

    output$plot2 <- renderPlot({
        req(consensusSelection())
        if (nrow(consensusSelection()) <= 51) {
            plotData(consensusSelection(), type = "nucleotide")
        } else {
            NULL
        }
    })

    output$plot3 <- renderPlot({
        req(oligoSelection())
        plotData(oligoSelection())
    })

    output$plot4 <- renderPlot({
        req(selectedOligo())
        from <- selectedOligo()$start
        to <- selectedOligo()$end
        plotData(consensus()[
            consensus()$position >= from & consensus()$position <= to,
            ], type = "nucleotide")
    })

    output$plot5 <- renderPlot({
        req(selectedOligoMatch())
        plotData(selectedOligoMatch())
    })

    #### Tables ####

    output$table1 <- DT::renderDataTable({
        req(consensusSelection())
        roundDbls(as.data.frame(consensusSelection()))
    }, options = list(
        info = FALSE,
        searching = FALSE, paging = FALSE,
        scrollX = TRUE, autoWidth = FALSE,
        ordering = FALSE, scrollY = "400"
        ), rownames = FALSE, selection  = "none"
    )

    output$table2 <- DT::renderDataTable({
        req(oligoSelection())
        roundDbls(removeListColumns(as.data.frame(oligoSelection())))
    }, options = list(
        info = FALSE,
        searching = FALSE, paging = FALSE,
        scrollX = TRUE, autoWidth = TRUE,
        ordering = TRUE, scrollY = "400"
       ), rownames = FALSE,
    selection = list(mode = "single", selected = 1)
    )

    output$table3 <- DT::renderDataTable({
        req(selectedOligo())
        roundDbls(as.data.frame(selectedOligo()))
    }, options = list(
        info = FALSE,
        searching = FALSE, paging = FALSE,
        scrollX = TRUE, autoWidth = TRUE,
        ordering = FALSE
    ), rownames = FALSE, selection  = "none"
    )

    output$table4 <-  DT::renderDataTable({
        req(selectedOligoMatch())
        x <- as.data.frame(selectedOligoMatch())
        x <- x[!grepl("id", names(x))]
        roundDbls(x)
    }, options = list(
        info = FALSE,
        searching = FALSE, paging = FALSE,
        scrollX = TRUE, autoWidth = FALSE,
        ordering = FALSE
    ), rownames = FALSE, selection  = "none"
    )

    #### Download links ####

    output$download1 <- downloadHandler(
        filename <- function() {
            paste0("consensus_profile-", Sys.Date(), ".csv")
        },
        content <- function(file) {
            write.csv(as.data.frame(consensusSelection()), file)
        })

    output$download2 <- downloadHandler(
        filename <- function() {
            paste0("oligos-", Sys.Date(), ".csv")
        },
        content <- function(file) {
            write.csv(as.data.frame(oligoSelection()), file)
        })


}

# Run app ======================================================================

shinyApp(ui, server)
