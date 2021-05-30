library(shiny)
library(shinydashboard)
library(shinycssloaders)
devtools::load_all("C:/Users/Sofia/Desktop/rprimer")

data("exampleRprimerAlignment")

plotWidth = 800
plotHeight = 600

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

# UI ===========================================================================

ui <- dashboardPage(
    dashboardHeader(title = ""),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Target alignment", tabName = "import"),
            menuItem("Consensus profile", tabName = "consensus"),
            menuItem("Oligos", tabName = "oligos"),
            menuItem("Assays", tabName = "assays")
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
                upload,
                useExample,
                hr(),
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
                column(width = 3,
                h4("Settings"),
                br(),
                numericInput(
                    "ambiguityThreshold", h5("Threshold for an ambiguous base"),
                    value = 0.05, min = 0, max = 0.2, width = 200
                ),
                br(),
                h5("Region of interest"),
                uiOutput("roiFrom"),
                uiOutput("roiTo"),
                br(),
                actionButton("getConsensusProfile", "Get consensus profile")
                ),
                column(width = 9,
                h4("Output"),
                br(),
                h5("Show region"),
                column(width = 6,
                uiOutput("zoomFrom")
                ),
                column(width = 6,
                uiOutput("zoomTo")
                ),
                br(),
                tabBox(title = "", width = 12,
                       tabPanel(title = "Plot",
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
            box(width = 12, title = "Oligos",
                h5("Design settings"),
                actionButton("getOligos", "Design oligos"),
                hr(),
                plotOutput("plot3", height = plotHeight, width = plotWidth)
            )
        ),

        #### Assays ####

        tabItem(tabName = "assays",
            box(width = 12, title = "Assays",
                h5("Design settings"),
                actionButton("getAssays", "Design assays"),
                hr(),
                plotOutput("plot4", height = plotHeight, width = plotWidth)
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
        oligos(consensus())
    })

    assayCandidates <- eventReactive(input$getAssays, {
        req(oligoCandidates())
        assays(oligoCandidates())
    })

    #### Render UI ####


    output$roiFrom <- renderUI({
        req(aln())
        numericInput(
            "roiFrom", h5("From"), min = 1, max = ncol(aln()) - 1, value = 1,
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


    #### Html ####

    output$html1 <- renderText({
        req(aln())
        paste("Number of sequences:", nrow(aln()))
    })

    output$html2 <- renderText({
        req(aln())
        paste("Alignment length:", ncol(aln()))
    })

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
        req(oligoCandidates())
        plotData(oligoCandidates())
    })

    output$plot4 <- renderPlot({
        req(assayCandidates())
        plotData(assayCandidates())
    })

    #### Tables ####

    output$table1 <- DT::renderDataTable({
        req(consensusSelection())
        as.data.frame(consensusSelection())
    }, options = list(
        searching = FALSE, paging = FALSE,
        scrollX = TRUE, autoWidth = TRUE,
        ordering = FALSE, scrollY = "600"), rownames = FALSE
    )

    #### Download links ####

    output$download1 <- downloadHandler(
        filename <- function() {
            paste0("consensus_profile-", Sys.Date(), ".csv")
        },
        content <- function(file) {
            write.csv(as.data.frame(consensusSelection()), file)
        })


}

# Run app ======================================================================

shinyApp(ui, server)
