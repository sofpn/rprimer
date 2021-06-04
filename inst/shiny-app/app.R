#library(shiny)
library(shinydashboard)
#library(shinycssloaders)
#devtools::load_all("C:/Users/Sofia/Desktop/rprimer")

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

makeListTable <- function(x) {
  x <- as.data.frame(x)
  list <- vapply(seq_len(ncol(x)), function(i) {
    is.list(x[, i])
  }, logical(1L))
  x <- x[, list]
  allCols <- lapply(seq_len(ncol(x)), function(i) {
    x[ , i][[1]]
  })
  all <- data.frame(do.call("cbind", allCols))
  names(all) <- names(x)
  all$gcContent <- as.double(all$gcContent)
  all$tm <- as.double(all$tm)
  all$deltaG <- as.double(all$deltaG)
  all
}

makeEmptyRow <- function(x) {
  emptyRow <- x[1, , drop = FALSE]
  emptyRow[1, ] <- NA
  emptyRow$roiStart <- x$roiStart[[1]]
  emptyRow$roiEnd <- x$roiEnd[[1]]
  emptyRow$start <- 1
  emptyRow$end <- 1
  emptyRow$fwd <- emptyRow$rev <- TRUE
  emptyRow$type <- "primer"
  emptyRow$iupacSequence <- "emptyRow"
  emptyRow
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
  h5("Hepatitis E virus alignment")
)

useProbe <- conditionalPanel(
  condition = "input.probe == true",
  sliderInput(
    "lengthProbe",
    h5("Length"),
    value = c(18, 22), min = 15, max = 40,

  ),
  numericInput(
    "maxDegeneracyProbe",
    h5("Maximum degeneracy"),
    value = 4, min = 1, max = 64,

  ),
  checkboxInput(
    "avoidFiveEndGProbe",
    h5("Avoid 5' end G"),
    value = TRUE,

  ),
  sliderInput(
    "gcProbe",
    h5("GC content range"),
    value = c(0.4, 0.65), min = 0, max = 1,

  ),
  sliderInput(
    "tmProbe",
    h5("Melting temperature range"),
    value = c(50, 70), min = 20, max = 90,

  ),
  numericInput(
    "concProbe",
    h5("Concentration (nM)"),
    value = 250, min = 20, max = 2000,

  )
)

# UI ===========================================================================

ui <- shinydashboard::dashboardPage(
  shinydashboard::dashboardHeader(title = ""),
  shinydashboard::dashboardSidebar(
    sidebarMenu(
      menuItem("Welcome!", tabName = "welcome"),
      menuItem("Target alignment", tabName = "import"),
      menuItem("1.  Consensus profile", tabName = "consensus"),
      menuItem("2.  Oligos", tabName = "oligos"),
      menuItem("3.  Assays", tabName = "assays")
    )
  ),
  dashboardBody(
    tabItems(

      #### Welcome ####

      tabItem(tabName = "welcome",
              fluidRow(
                box(
                  title = ("Instructions for use"),
                  h5(""),
                  br(),
                  uiOutput("codelink"),
                  width = 12
                )
              )
      ),

      #### Target alignment ####

      tabItem(tabName = "import",
              fluidRow(
                column(width = 12,
                       box(width = 12, title = "Select target alignment",
                           radioButtons(
                             "fileInput", label = NULL,
                             choices = c(
                               "Upload file", "Use example data"
                             ),
                             selected = "Upload file"),
                           upload,
                           useExample,
                           br(),
                           htmlOutput("html1"),
                           htmlOutput("html2")
                       )
                ))
      ),

      #### Consensus ####

      tabItem(tabName = "consensus",
              fluidRow(
                column(width = 3,
                       box(width = 12, title = "Settings",
                           numericInput(
                             "ambiguityThreshold",
                             h5("Threshold for an ambiguous base"),
                             value = 0.05, min = 0, max = 0.2,
                           ),
                           uiOutput("roi"),
                           uiOutput("roiFrom"),
                           uiOutput("roiTo"),
                           br(),
                           uiOutput("getConsensusProfile")
                       )),
                column(width =  9,
                       box(width = 12, title = "Output",
                           uiOutput("showRegion"),
                           column(width = 2,
                                  uiOutput("zoomFrom")
                           ),
                           column(width = 2,
                                  uiOutput("zoomTo")
                           ),
                           br(),
                           tabBox(title = "",
                                  width = 12,
                                  tabPanel(title = "Plot",
                                           br(),
                                           shinycssloaders::withSpinner(
                                             plotOutput(
                                               "plot1",
                                               height = plotHeight,
                                               width = "100%"
                                             ), color = "grey")
                                  ),
                                  tabPanel(title = "Table",
                                           br(),
                                           downloadLink(
                                             "download1", "Download"
                                           ),
                                           br(),
                                           br(),
                                           DT::dataTableOutput("table1")
                                  )
                           ))))),

      #### Oligos ####

      tabItem(tabName = "oligos",
              fluidRow(column(width = 3,
                              tags$head(tags$style(
                                HTML(
                                  '.box{-webkit-box-shadow: none;
                    -moz-box-shadow: none;
                    box-shadow: none;}'))),
                              box(width = 12, title = "Settings",
                                  h4("Primers"),
                                  hr(),
                                  sliderInput(
                                    "lengthPrimer",
                                    h5("Length"),
                                    value = c(18, 22), min = 15, max = 40,

                                  ),
                                  numericInput(
                                    "maxDegeneracyPrimer",
                                    h5("Maximum degeneracy"),
                                    value = 4, min = 1, max = 64,

                                  ),
                                  checkboxInput(
                                    "avoidThreeEndRunsPrimer",
                                    h5("Avoid 3' end runs"),
                                    value = TRUE,

                                  ),
                                  checkboxInput(
                                    "gcClampPrimer",
                                    h5("Use GC clamp"),
                                    value = TRUE,

                                  ),
                                  sliderInput(
                                    "gcPrimer",
                                    h5("GC content range"),
                                    value = c(0.4, 0.65), min = 0, max = 1,

                                  ),
                                  sliderInput(
                                    "tmPrimer",
                                    h5("Melting temperature range"),
                                    value = c(50, 65), min = 20, max = 90,

                                  ),
                                  numericInput(
                                    "concPrimer",
                                    h5("Concentration (nM)"),
                                    value = 500, min = 20, max = 2000,

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

                                  ),
                                  useProbe,
                                  br(),
                                  h4("General"),
                                  hr(),
                                  numericInput(
                                    "maxGapFrequency",
                                    h5("Maximum gap frequency"),
                                    value = 0.01, min = 0, max = 0.2,

                                  ),
                                  numericInput(
                                    "concNa",
                                    h5("Sodium ion concentration (M)"),
                                    value = 0.05, min = 0, max = 1,

                                  ),
                                  hr(),
                                  uiOutput("getOligos")
                              )),
                       column(width = 9,
                              box(width = 12, title = "Filter",
                                  column(width = 2, uiOutput("fwdRegionFrom")),
                                  column(width = 2, uiOutput("fwdRegionTo")),
                                  column(width = 2, uiOutput("revRegionFrom")),
                                  column(width = 2, uiOutput("revRegionTo")),
                                  column(width = 2, uiOutput("prRegionFrom")),
                                  column(width = 2, uiOutput("prRegionTo"))
                              ),
                              box(width = 12, title = "All oligos",
                                  tabBox(title = "", width = 12,
                                         tabPanel(title = "Plot",
                                                  br(),
                                                  shinycssloaders::withSpinner(
                                                    plotOutput(
                                                      "plot3",
                                                      height = plotHeight,
                                                      width = "100%"
                                                    ), color = "grey")
                                         ),
                                         tabPanel(title = "Table and selection",
                                                  downloadLink(
                                                    "download2", "Download table"
                                                  ),
                                                  br(),
                                                  downloadLink(
                                                    "download3", "Download fasta"
                                                  ),
                                                  br(),
                                                  br(),
                                                  DT::dataTableOutput("table2")
                                         ))),
                              box(width = 12, title = "Selection",
                                  collapsible = TRUE, collapsed = FALSE,
                                  box(width = 12, title = "Table",
                                      solidHeader = TRUE,
                                      downloadLink(
                                        "download4", "Download table"
                                      ),
                                      br(),
                                      br(),
                                      DT::dataTableOutput("table3")
                                  ),
                                  box(width = 12, title = "All sequence variants",
                                      solidHeader = TRUE,
                                      downloadLink(
                                        "download5", "Download fasta"
                                      ),
                                      br(),
                                      br(),
                                      DT::dataTableOutput("table4")
                                  ),
                                  box(width = 6,
                                      title = "Nucleotide distribution in target alignment",
                                      solidHeader = TRUE,
                                      shinycssloaders::withSpinner(plotOutput(
                                        "plot4",
                                        height = plotHeight / 2,
                                        width = "100%"
                                      ), color = "grey")
                                  ),
                                  box(width = 6, title = "Match plot",
                                      solidHeader = TRUE,
                                      shinycssloaders::withSpinner(plotOutput(
                                        "plot5",
                                        height = plotHeight / 2,
                                        width = "100%"
                                      ), color = "grey")
                                  ),
                                  box(width = 12, title = "Match details",
                                      solidHeader = TRUE,
                                      htmlOutput("html3"),
                                      br(),
                                      br(),
                                      htmlOutput("html4"),
                                      br(),
                                      br(),
                                      htmlOutput("html5"),
                                      br(),
                                      br(),
                                      htmlOutput("html6"),
                                      br(),
                                      br(),
                                      htmlOutput("html7"),
                                      br(),
                                      br(),
                                      htmlOutput("html8")
                                  ))



                       ))),

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
           gcPrimer = input$gcPrimer,
           tmPrimer = input$tmPrimer,
           concPrimer = input$concPrimer,
           designStrategyPrimer = input$designStrategyPrimer,
           probe = input$probe,
           lengthProbe = input$lengthProbe,
           maxDegeneracyProbe = input$maxDegeneracyProbe,
           avoidFiveEndGProbe = input$avoidFiveEndGProbe,
           gcProbe = input$gcProbe,
           tmProbe = input$tmProbe,
           concNa = input$concNa
    )
  })

  oligoSelection <- reactive({
    req(oligoCandidates())
    x <- oligoCandidates()
    x <- as.data.frame(x)
    emptyRow <- makeEmptyRow(x)
    fwd <- x$fwd & x$type == "primer" &
      x$start >= as.numeric(input$fwdRegionFrom)  &
      x$end <= as.numeric(input$fwdRegionTo)
    rev <- (x$rev & x$type == "primer" &
              x$start >= as.numeric(input$revRegionFrom)  &
              x$end <= as.numeric(input$revRegionTo))
    if (any(x$type == "probe")) {
      pr <- x$type == "probe" &
        x$start >= as.numeric(input$prRegionFrom)  &
        x$end <= as.numeric(input$prRegionTo)
      x <- x[fwd | rev | pr, ]
    } else {
      x <- x[fwd | rev, ]
    }
    x$fwd[x$type == "primer" & x$start < as.numeric(input$fwdRegionFrom)] <- FALSE
    x$fwd[x$type == "Primer" & x$end > as.numeric(input$fwdRegionTo)] <- FALSE
    x$rev[x$type == "primer" & x$start < as.numeric(input$revRegionFrom)] <- FALSE
    x$rev[x$type == "primer" & x$end > as.numeric(input$revRegionTo)] <- FALSE
    if (nrow(x) == 0L) {
      x <- emptyRow
    }
    RprimerOligo(x)
  })

  selectedOligo <- reactive({
    req(oligoSelection())
    if (!is.null(input$table2_rows_selected)) {
      oligoSelection()[input$table2_rows_selected, ]
    } else {
      oligoSelection()[1, ]
    }
  })

  selectedOligoMatch <- reactive({
    req(selectedOligo())
    req(aln())
    if (selectedOligo()$iupacSequence == "emptyRow") {
      NULL
    } else {
      checkMatch(selectedOligo(), aln())
    }
  })

  assayCandidates <- eventReactive(input$getAssays, {
    req(oligoCandidates())
    assays(oligoCandidates())
  })

  #### Render UI ####

  url <- a("Source code", href = "https://github.com/sofpn/rprimer")
  output$codelink <- renderUI({
    tagList(url)
  })

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

    )
  })

  output$roiTo <- renderUI({
    req(aln())
    numericInput(
      "roiTo", h5("To"), min = 2, max = ncol(aln()), value = ncol(aln()),

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
    if (any(oligoCandidates()$type == "probe")) {
      numericInput("prRegionFrom", h5("Probe, from"),
                   min = oligoCandidates()$roiStart[[1]],
                   max = oligoCandidates()$roiEnd[[1]],
                   value = oligoCandidates()$roiStart[[1]]
      )
    } else {
      NULL
    }

  })

  output$prRegionTo <- renderUI({
    req(oligoCandidates())
    if (any(oligoCandidates()$type == "probe")) {
      numericInput("prRegionTo", h5("Probe, to"),
                   min = oligoCandidates()$roiStart[[1]],
                   max = oligoCandidates()$roiEnd[[1]],
                   value = oligoCandidates()$roiEnd[[1]]
      )
    } else {
      NULL
    }

  })

  output$allVariants <- renderUI({
    req(selectedOligo())
    h4("All sequence variants")
  })


  output$ntDistribution <- renderUI({
    req(selectedOligo())
    h4("Nucleotide distribution in target alignment")
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

  output$html3 <- renderText({
    req(!is.null(selectedOligoMatch()))
    c(
      "<b>Perfect match</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedOligoMatch()$perfectMatch, 2), "<br><br>",
      "ID<br>",
      selectedOligoMatch()$idPerfectMatch[[1]])
  })

  output$html4 <- renderText({
    req(!is.null(selectedOligoMatch()))
    c(
      "<b>One mismatch</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedOligoMatch()$oneMismatch, 2), "<br><br>",
      "ID<br>",
      selectedOligoMatch()$idOneMismatch[[1]])
  })

  output$html5 <- renderText({
    req(!is.null(selectedOligoMatch()))
    c(
      "<b>Two mismatches</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedOligoMatch()$twoMismatches, 2), "<br><br>",
      "ID<br>",
      selectedOligoMatch()$idTwoMismatches[[1]])
  })

  output$html6 <- renderText({
    req(!is.null(selectedOligoMatch()))
    c(
      "<b>Three mismatches</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedOligoMatch()$threeMismatches, 2), "<br><br>",
      "ID<br>",
      selectedOligoMatch()$idThreeMismatches[[1]])
  })

  output$html7 <- renderText({
    req(!is.null(selectedOligoMatch()))
    c(
      "<b>Four or more mismatches</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedOligoMatch()$fourOrMoreMismatches, 2), "<br><br>",
      "ID<br>",
      selectedOligoMatch()$idFourOrMoreMismatches[[1]])
  })

  output$html8 <- renderText({
    req(!is.null(selectedOligoMatch()))
    c(
      "<b>Off target match</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedOligoMatch()$offTargetMatch, 2), "<br><br>",
      "ID<br>",
      selectedOligoMatch()$idOffTargetMatch[[1]])
  })




  #### Plots ####

  output$plot1 <- renderPlot({
    req(consensusSelection())
    plotData(consensusSelection())
  })

  output$plot3 <- renderPlot({
    req(oligoSelection())
    plotData(oligoSelection())
  })

  output$plot4 <- renderPlot({
    req(selectedOligo())
    if (selectedOligo()$iupacSequence == "emptyRow") {
      NULL
    } else {
      from <- selectedOligo()$start
      to <- selectedOligo()$end
      plotData(consensus()[
        consensus()$position >= from & consensus()$position <= to,
      ], type = "nucleotide")
    }
  })

  output$plot5 <- renderPlot({
    req(!is.null(selectedOligoMatch()))
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
    ordering = FALSE, scrollY = "300"
  ), rownames = FALSE, selection  = "none"
  )

  output$table2 <- DT::renderDataTable({
    req(oligoSelection())
    if (oligoSelection()$iupacSequence[[1]] == "emptyRow") {
      NULL
    } else {
      roundDbls(removeListColumns(as.data.frame(oligoSelection())))
    }
  }, options = list(
    info = FALSE,
    searching = FALSE, paging = FALSE,
    scrollX = TRUE, autoWidth = TRUE,
    ordering = TRUE, scrollY = "300"
  ), rownames = FALSE,
  selection = list(mode = "single", selected = 1)
  )

  output$table3 <- DT::renderDataTable({
    req(selectedOligo())
    if (selectedOligo()$iupacSequence[[1]] == "emptyRow") {
      NULL
    } else {
      x <- selectedOligo()
      x <- as.data.frame(x)
      x <- removeListColumns(x)
      roundDbls(x)
    }
  }, options = list(
    info = FALSE,
    searching = FALSE, paging = FALSE,
    scrollX = TRUE, autoWidth = TRUE,
    ordering = FALSE
  ), rownames = FALSE, selection  = "none"
  )

  output$table4 <- DT::renderDataTable({
    req(selectedOligo())
    if (selectedOligo()$iupacSequence[[1]] == "emptyRow") {
      NULL
    } else {
      roundDbls(makeListTable(selectedOligo()))
    }
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

  output$download3 <- downloadHandler(
    filename <- function() {
      paste0("oligos-fasta", Sys.Date(), ".txt")
    },
    content <- function(file) {
      x <- oligoSelection()
      x <- as(oligoSelection(), "DNAStringSet")
      Biostrings::writeXStringSet(x, file)
    })

  output$download4 <- downloadHandler(
    filename <- function() {
      paste0("oligo-selection-", Sys.Date(), ".csv")
    },
    content <- function(file) {
      write.csv(as.data.frame(selectedOligo()), file)
    })

  output$download5 <- downloadHandler(
    filename <- function() {
      paste0("oligo-selection-fasta-", Sys.Date(), ".txt")
    },
    content <- function(file) {
      x <- selectedOligo()
      x <- as(selectedOligo(), "DNAStringSet")
      Biostrings::writeXStringSet(x, file)
    })
}

# Run app ======================================================================

shinyApp(ui, server)
