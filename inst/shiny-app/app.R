#library(shiny)
#library(shinydashboard)
#devtools::load_all("C:/Users/Sofia/Desktop/rprimer")

# remove ws  plots o match info pr
# make selection renderUI, or sep boxes

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


splitAssayToList <- function(x) {
  x <- as.data.frame(x)
  fwd <- x[, grepl("Fwd", names(x))]
  names(fwd) <- gsub("Fwd", "", names(fwd))
  rev <- x[, grepl("Rev", names(x))]
  names(rev) <- gsub("Rev", "", names(rev))
  if (any(grepl("Pr", names(x)))) {
    pr <- x[, grepl("Pr", names(x))]
    names(pr) <- gsub("Pr", "", names(pr))
    all <- list(fwd, rev, pr)
  } else {
    all <- list(fwd, rev)
  }
  all
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
  emptyRow$length <- NA
  emptyRow
}

# Conditional panels ===========================================================

upload <- conditionalPanel(
  condition = "input.fileInput == 'Upload file'",
  br(),
  h5("Upload a multiple DNA sequence alignment"),
  fileInput(
    "file1", "", multiple = FALSE, accept = c("text")
  ),
  radioButtons(
    "filetype", h5("Format"),
    choices = c("Fasta" = "fasta", "Clustal" = "clustal"),
    selected = "fasta"
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
    h5("Maximum degeneracy (1-64)"),
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
    h5("Melting temperature range (Celcius degrees)"),
    value = c(50, 70), min = 20, max = 90,

  ),
  numericInput(
    "concProbe",
    h5("Concentration (nM)"),
    value = 250, min = 20, max = 2000,

  )
)

# UI ===========================================================================

ui <- dashboardPage(
  dashboardHeader(title = "rprimer"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Welcome", tabName = "welcome", icon = icon("heart")),
      menuItem("Target alignment", tabName = "import", icon = icon("file-upload")),
      menuItem("Consensus profile", tabName = "consensus", icon = icon("th")),
      menuItem("Oligos", tabName = "oligos", icon = icon("dna")),
      menuItem("Assays", tabName = "assays", icon = icon("stream"))
    )
  ),
  dashboardBody(
    tabItems(

      #### Welcome ####

      tabItem(tabName = "welcome",
              fluidRow(
                box(
                  title = ("Design degenerate oligos from a multiple DNA sequence alignment"),
                  h5(tags$b("Introduction")),
                  h5(
                  "rprimer provides tools for visualizing sequence
                  conservation and designing degenerate primers, probes and (RT)-(q/d)PCR
                  assays from a multiple DNA sequence alignment. The workflow is
                  developed primarily for sequence variable RNA viruses, but it should be
                  equally useful for other organisms with high sequence variability"
                  ),
                  br(),
                  h5(tags$b("Instructions for use")),
                  h5("Add link to package vignette"),
                  br(),
                  h5(tags$b("Citation")),
                  h5("S Persson et al., manuscript in preparation"),
                  br(),
                  h5(tags$b("Source code")),
                  uiOutput("codelink"),
                  br(),
                  h5(tags$b("Contact")),
                  h5("sofia.persson@slv.se"),
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
                                    h5("Maximum degeneracy (1-64)"),
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
                                    h5("Melting temperature range (Celcius degrees)"),
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
                                  box(width = 12,
                                      title = "Nucleotide distribution in target alignment",
                                      solidHeader = TRUE,
                                      shinycssloaders::withSpinner(plotOutput(
                                        "plot4",
                                        height = plotHeight / 2,
                                        width = "100%"
                                      ), color = "grey")
                                  ),
                                  box(width = 12, title = "Match plot",
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
                                      htmlOutput("html7")
                                  ))



                       ))),

      #### Assays ####

      tabItem(tabName = "assays",
              fluidRow(column(width = 3,
                              tags$head(tags$style(
                                HTML(
                                  '.box{-webkit-box-shadow: none;
                    -moz-box-shadow: none;
                    box-shadow: none;}'))),
                              box(width = 12, title = "Settings",
                                  sliderInput(
                                    "length",
                                    h5("Length"),
                                    value = c(60, 120), min = 40, max = 5000,
                                  ),
                                  numericInput(
                                    "tmDifferencePrimers",
                                    h5("Maximum melting temperature difference between primers (Celcius degrees)"),
                                    value = 10, min = 0, max = Inf,
                                  ),
                                  hr(),
                                  uiOutput("getAssays")
                              )),
                       column(width = 9,
                              box(width = 12, title = "Filter",
                                  column(width = 2, uiOutput("assayRegionFrom")),
                                  column(width = 2, uiOutput("assayRegionTo")),
                              ),
                              box(width = 12, title = "All assays",
                                  tabBox(title = "", width = 12,
                                         tabPanel(title = "Plot",
                                                  br(),
                                                  shinycssloaders::withSpinner(
                                                    plotOutput(
                                                      "plot6",
                                                      height = plotHeight,
                                                      width = "100%"
                                                    ), color = "grey")
                                         ),
                                         tabPanel(title = "Table and selection",
                                                  downloadLink(
                                                    "download6", "Download table"
                                                  ),
                                                  br(),
                                                  downloadLink(
                                                    "download7", "Download fasta"
                                                  ),
                                                  br(),
                                                  br(),
                                                  DT::dataTableOutput("table5")
                                         ))),
                              box(width = 12, title = "Selection",
                                  collapsible = TRUE, collapsed = FALSE,
                                  box(width = 12, title = "Table",
                                      solidHeader = TRUE,
                                      downloadLink(
                                        "download8", "Download table"
                                      ),
                                      br(),
                                      br(),
                                      DT::dataTableOutput("table6")
                                  ),
                                  box(width = 12, title = "All sequence variants",
                                      solidHeader = TRUE,
                                      downloadLink(
                                        "download9", "Download fasta"
                                      ),
                                      br(),
                                      br(),
                                      uiOutput("fwd2"),
                                      DT::dataTableOutput("table7"),
                                      br(),
                                      uiOutput("rev2"),
                                      DT::dataTableOutput("table8"),
                                      br(),
                                      uiOutput("pr2"),
                                      DT::dataTableOutput("table9")

                                  ),
                                  box(width = 12,
                                      title = "Nucleotide distribution in target alignment",
                                      solidHeader = TRUE,
                                      uiOutput("fwd"),
                                      br(),
                                      br(),
                                      shinycssloaders::withSpinner(plotOutput(
                                        "plot7",
                                        height = plotHeight / 2,
                                        width = "100%"
                                      ), color = "grey"),
                                      uiOutput("rev"),
                                      br(),
                                      br(),
                                      shinycssloaders::withSpinner(plotOutput(
                                        "plot8",
                                        height = plotHeight / 2,
                                        width = "100%"
                                      ), color = "grey"),
                                      uiOutput("pr"),
                                      br(),
                                      br(),
                                      shinycssloaders::withSpinner(plotOutput(
                                        "plot9",
                                        height = plotHeight / 2,
                                        width = "100%"
                                      ), color = "grey")
                                  ),
                                  box(width = 12, title = "Match plot",
                                      solidHeader = TRUE,
                                      shinycssloaders::withSpinner(plotOutput(
                                        "plot10",
                                        height = plotHeight / 2,
                                        width = "100%"
                                      ), color = "grey")
                                  ),
                                  box(width = 12, title = "Match details",
                                      solidHeader = TRUE,
                                      uiOutput("fwd3"),
                                      hr(),
                                      htmlOutput("html9"),
                                      br(),
                                      br(),
                                      htmlOutput("html10"),
                                      br(),
                                      br(),
                                      htmlOutput("html11"),
                                      br(),
                                      br(),
                                      htmlOutput("html12"),
                                      br(),
                                      br(),
                                      htmlOutput("html13"),
                                      br(),
                                      br(),
                                      uiOutput("rev3"),
                                      hr(),
                                      htmlOutput("html14"),
                                      br(),
                                      br(),
                                      htmlOutput("html15"),
                                      br(),
                                      br(),
                                      htmlOutput("html16"),
                                      br(),
                                      br(),
                                      htmlOutput("html17"),
                                      br(),
                                      br(),
                                      htmlOutput("html18"),
                                      br(),
                                      br(),
                                      uiOutput("pr3"),
                                      hr(),
                                      htmlOutput("html19"),
                                      br(),
                                      br(),
                                      htmlOutput("html20"),
                                      br(),
                                      br(),
                                      htmlOutput("html21"),
                                      br(),
                                      br(),
                                      htmlOutput("html22"),
                                      br(),
                                      br(),
                                      htmlOutput("html23")
                                  ))


                       )



              )))))

# Server =======================================================================

server <- function(input, output) {

  #### Data ####

  aln <- reactive({
    if (input$fileInput == "Upload file") {
      req(input$file1)
      x <- Biostrings::readDNAMultipleAlignment(
        input$file1$datapath, format = input$filetype
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
    if (is.na(selectedOligo()$length[[1]])) {
      NULL
    } else {
      checkMatch(selectedOligo(), aln())
    }
  })

  assayCandidates <- eventReactive(input$getAssays, {
    req(oligoCandidates())
    assays(oligoCandidates(),
           length = input$length,
           tmDifferencePrimers = input$tmDifferencePrimers)
  })

  assaySelection <- reactive({
    req(assayCandidates())
    x <- assayCandidates()
    x <- as.data.frame(x)
    emptyRow <- makeEmptyRow(x)
    x <- x[
      x$start >= input$assayRegionFrom & x$end <= input$assayRegionTo,
    ]
    if (nrow(x) == 0L) {
      x <- emptyRow
    }
    RprimerAssay(x)
  })

  selectedAssay <- reactive({
    req(assaySelection())
    if (!is.null(input$table5_rows_selected)) {
      assaySelection()[input$table5_rows_selected, ]
    } else {
      assaySelection()[1, ]
    }
  })

  selectedAssayList <- reactive({
    req(!is.na(selectedAssay()$length[[1]]))
    splitAssayToList(selectedAssay())
  })

  selectedAssayMatch <- reactive({
    req(selectedAssay())
    req(aln())
    if (is.na(selectedAssay()$length[[1]])) {
      NULL
    } else {
      checkMatch(selectedAssay(), aln())
    }
  })

  #### Render UI ####

  url <- a("GitHub", href = "https://github.com/sofpn/rprimer")
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

  output$getAssays <- renderUI({
    req(oligoSelection())
    actionButton("getAssays", "Get assays")
  })

  output$assayRegionFrom <- renderUI({
    req(assayCandidates())
    numericInput("assayRegionFrom", h5("From"),
                 min = assayCandidates()$roiStart[[1]],
                 max = assayCandidates()$roiEnd[[1]],
                 value = assayCandidates()$roiStart[[1]]
    )
  })

  output$assayRegionTo <- renderUI({
    req(assayCandidates())
    numericInput("assayRegionTo", h5("To"),
                 min = assayCandidates()$roiStart[[1]],
                 max = assayCandidates()$roiEnd[[1]],
                 value = assayCandidates()$roiEnd[[1]]
    )
  })

  output$fwd <- renderUI({
    req(selectedAssay())
    h4("Forward")
  })

  output$rev <- renderUI({
    req(selectedAssay())
    h4("Reverse")
  })

  output$pr <- renderUI({
    req(any(grepl("Pr", names(selectedAssay()))))
    h4("Probe")
  })

  output$fwd2 <- renderUI({
    req(selectedAssay())
    h4("Forward")
  })

  output$rev2 <- renderUI({
    req(selectedAssay())
    h4("Reverse")
  })

  output$pr2 <- renderUI({
    req(any(grepl("Pr", names(selectedAssay()))))
    h4("Probe")
  })

  output$fwd3 <- renderUI({
    req(selectedAssayMatch())
    h4("Forward")
  })

  output$rev3 <- renderUI({
    req(selectedAssayMatch())
    h4("Reverse")
  })

  output$pr3 <- renderUI({
    req(any(grepl("Pr", names(selectedAssayMatch()))))
    h4("Probe")
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

  output$html9 <- renderText({
    req(!is.null(selectedAssayMatch()))
    c(
      "<b>Perfect match</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedAssayMatch()$perfectMatchFwd, 2), "<br><br>",
      "ID<br>",
      selectedAssayMatch()$idPerfectMatchFwd[[1]])
  })

  output$html10 <- renderText({
    req(!is.null(selectedAssayMatch()))
    c(
      "<b>One mismatch</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedAssayMatch()$oneMismatchFwd, 2), "<br><br>",
      "ID<br>",
      selectedAssayMatch()$idOneMismatchFwd[[1]])
  })

  output$html11 <- renderText({
    req(!is.null(selectedAssayMatch()))
    c(
      "<b>Two mismatches</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedAssayMatch()$twoMismatchesFwd, 2), "<br><br>",
      "ID<br>",
      selectedAssayMatch()$idTwoMismatchesFwd[[1]])
  })

  output$html12 <- renderText({
    req(!is.null(selectedAssayMatch()))
    c(
      "<b>Three mismatches</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedAssayMatch()$threeMismatchesFwd, 2), "<br><br>",
      "ID<br>",
      selectedAssayMatch()$idThreeMismatchesFwd[[1]])
  })

  output$html13 <- renderText({
    req(!is.null(selectedAssayMatch()))
    c(
      "<b>Four or more mismatches</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedAssayMatch()$fourOrMoreMismatchesFwd, 2), "<br><br>",
      "ID<br>",
      selectedAssayMatch()$idFourOrMoreMismatchesFwd[[1]])
  })

  output$html14 <- renderText({
    req(!is.null(selectedAssayMatch()))
    c(
      "<b>Perfect match</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedAssayMatch()$perfectMatchRev, 2), "<br><br>",
      "ID<br>",
      selectedAssayMatch()$idPerfectMatchRev[[1]])
  })

  output$html15 <- renderText({
    req(!is.null(selectedAssayMatch()))
    c(
      "<b>One mismatch</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedAssayMatch()$oneMismatchRev, 2), "<br><br>",
      "ID<br>",
      selectedAssayMatch()$idOneMismatchRev[[1]])
  })

  output$html16 <- renderText({
    req(!is.null(selectedAssayMatch()))
    c(
      "<b>Two mismatches</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedAssayMatch()$twoMismatchesRev, 2), "<br><br>",
      "ID<br>",
      selectedAssayMatch()$idTwoMismatchesRev[[1]])
  })

  output$html17 <- renderText({
    req(!is.null(selectedAssayMatch()))
    c(
      "<b>Three mismatches</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedAssayMatch()$threeMismatchesRev, 2), "<br><br>",
      "ID<br>",
      selectedAssayMatch()$idThreeMismatchesRev[[1]])
  })

  output$html18 <- renderText({
    req(!is.null(selectedAssayMatch()))
    c(
      "<b>Four or more mismatches</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedAssayMatch()$fourOrMoreMismatchesRev, 2), "<br><br>",
      "ID<br>",
      selectedAssayMatch()$idFourOrMoreMismatchesRev[[1]])
  })

  output$html19 <- renderText({
    req(any(grepl("Pr", names(selectedAssayMatch()))))
    c(
      "<b>Perfect match</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedAssayMatch()$perfectMatchPr, 2), "<br><br>",
      "ID<br>",
      selectedAssayMatch()$idPerfectMatchPr[[1]])
  })

  output$html20 <- renderText({
    req(any(grepl("Pr", names(selectedAssayMatch()))))
    c(
      "<b>One mismatch</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedAssayMatch()$oneMismatchPr, 2), "<br><br>",
      "ID<br>",
      selectedAssayMatch()$idOneMismatchPr[[1]])
  })

  output$html21 <- renderText({
    req(any(grepl("Pr", names(selectedAssayMatch()))))
    c(
      "<b>Two mismatches</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedAssayMatch()$twoMismatchesPr, 2), "<br><br>",
      "ID<br>",
      selectedAssayMatch()$idTwoMismatchesPr[[1]])
  })

  output$html22 <- renderText({
    req(any(grepl("Pr", names(selectedAssayMatch()))))
    c(
      "<b>Three mismatches</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedAssayMatch()$threeMismatchesPr, 2), "<br><br>",
      "ID<br>",
      selectedAssayMatch()$idThreeMismatchesPr[[1]])
  })

  output$html23 <- renderText({
    req(any(grepl("Pr", names(selectedAssayMatch()))))
    c(
      "<b>Four or more mismatches</b><br><br>",
      "Proportion of target sequences<br>",
      round(selectedAssayMatch()$fourOrMoreMismatchesPr, 2), "<br><br>",
      "ID<br>",
      selectedAssayMatch()$idFourOrMoreMismatchesPr[[1]])
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
    if (is.na(selectedOligo()$length[[1]])) {
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

  output$plot6 <- renderPlot({
    req(assaySelection())
    plotData(assaySelection())
  })

  output$plot7 <- renderPlot({
    req(selectedAssay())
    if (is.na(selectedAssay()$length[[1]])) {
      NULL
    } else {
      from <- selectedAssay()$startFwd
      to <- selectedAssay()$endFwd
      plotData(consensus()[
        consensus()$position >= from & consensus()$position <= to,
      ], type = "nucleotide")
    }
  })

  output$plot8 <- renderPlot({
    req(selectedAssay())
    if (is.na(selectedAssay()$length[[1]])) {
      NULL
    } else {
      from <- selectedAssay()$startRev
      to <- selectedAssay()$endRev
      plotData(consensus()[
        consensus()$position >= from & consensus()$position <= to,
      ], type = "nucleotide", rc = TRUE)
    }
  })

  output$plot9 <- renderPlot({
    req(any(grepl("Pr", names(selectedAssay()))))
    if (is.na(selectedAssay()$length[[1]])) {
      NULL
    } else {
      from <- selectedAssay()$startPr
      to <- selectedAssay()$endPr
      plotData(consensus()[
        consensus()$position >= from & consensus()$position <= to,
      ], type = "nucleotide", rc = TRUE)
    }
  })

  output$plot10 <- renderPlot({
    req(!is.null(selectedAssayMatch()))
    plotData(selectedAssayMatch())
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
    if (is.na(oligoSelection()$length[[1]])) {
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
    if (is.na(selectedOligo()$length[[1]])) {
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
    if (is.na(selectedOligo()$length[[1]])) {
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

  output$table5 <- DT::renderDataTable({
    req(assaySelection())
    if (is.na(assaySelection()$length[[1]])) {
      NULL
    } else {
      roundDbls(removeListColumns(as.data.frame(assaySelection())))
    }
  }, options = list(
    info = FALSE,
    searching = FALSE, paging = FALSE,
    scrollX = TRUE, autoWidth = TRUE,
    ordering = TRUE, scrollY = "300"
  ), rownames = FALSE,
  selection = list(mode = "single", selected = 1)
  )

  output$table6 <- DT::renderDataTable({
    req(selectedAssay())
    if (is.na(selectedAssay()$length[[1]])) {
      NULL
    } else {
      x <- selectedAssay()
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

  output$table7 <- DT::renderDataTable({
    req(selectedAssayList())
    roundDbls(makeListTable(selectedAssayList()[[1]]))
  }, options = list(
    info = FALSE,
    searching = FALSE, paging = FALSE,
    scrollX = TRUE, autoWidth = FALSE,
    ordering = FALSE
  ), rownames = FALSE, selection  = "none"
  )

  output$table8 <- DT::renderDataTable({
    req(selectedAssayList())
    roundDbls(makeListTable(selectedAssayList()[[2]]))
  }, options = list(
    info = FALSE,
    searching = FALSE, paging = FALSE,
    scrollX = TRUE, autoWidth = FALSE,
    ordering = FALSE
  ), rownames = FALSE, selection  = "none"
  )

  output$table9 <- DT::renderDataTable({
    req(length(selectedAssayList()) == 3)
    roundDbls(makeListTable(selectedAssayList()[[3]]))
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


  output$download6 <- downloadHandler(
    filename <- function() {
      paste0("assays-", Sys.Date(), ".csv")
    },
    content <- function(file) {
      write.csv(as.data.frame(assaySelection()), file)
    })

  output$download7 <- downloadHandler(
    filename <- function() {
      paste0("assays-fasta", Sys.Date(), ".txt")
    },
    content <- function(file) {
      x <- assayection()
      x <- as(assaySelection(), "DNAStringSet")
      Biostrings::writeXStringSet(x, file)
    })

  output$download8 <- downloadHandler(
    filename <- function() {
      paste0("assay-selection-", Sys.Date(), ".csv")
    },
    content <- function(file) {
      write.csv(as.data.frame(selectedAssay()), file)
    })

  output$download9 <- downloadHandler(
    filename <- function() {
      paste0("assay-selection-fasta-", Sys.Date(), ".txt")
    },
    content <- function(file) {
      x <- selectedAssay()
      x <- as(selectedAssay(), "DNAStringSet")
      Biostrings::writeXStringSet(x, file)
    })

}

# Run app ======================================================================

shinyApp(ui, server)
