library(shiny)
library(shinydashboard)
library(rprimer)
source("custom_functions.R")

data("exampleRprimerAlignment")

plotWidth = 800
plotHeight = 600

options(shiny.sanitize.errors = FALSE)

# Conditional panels ===========================================================

#### File upload/use example #####

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

#### Use probe ####

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

#### Assay with probe ####

probeSelectionMatch <- conditionalPanel(
    condition = "input.probe == true",
    br(),
    br(),
    h5("Probe"),
    br(),
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
)

assayProbeInfo <- conditionalPanel(
    condition = "input.probe == true",
    box(width = 12, title = "Probe",
        solidHeader = TRUE,
        DT::dataTableOutput("tablePr"),
        br(),
        h5("All sequence variants"),
        DT::dataTableOutput("table9"),
        br(),
        h5("Match information"),
        DT::dataTableOutput("tablePrMatch"),
        br(),
        h5("Nucleotide distribution in target alignment"),
        shinycssloaders::withSpinner(plotOutput(
            "plot9",
            height = plotHeight / 2,
            width = "100%"
        ), color = "grey"),
    )

)

# UI ===========================================================================

# Tab titles ===================================================================

ui <- dashboardPage(
    dashboardHeader(title = ""),
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

            # Welcome ==========================================================

            tabItem(tabName = "welcome",
                    fluidRow(
                        tags$head(tags$style(
                            HTML(
                                '.box{-webkit-box-shadow: none;
                    -moz-box-shadow: none;
                    box-shadow: none;}'))),
                        box(width = 12,
                            title = ("Design degenerate oligos from a multiple DNA sequence alignment"),
                            h5(
                                em("This application is under development")
                                ),
                            br(),
                            h5(tags$b("Introduction")),
                            h5(
                                "This application provides tools for visualizing sequence
                                conservation and designing degenerate primers, probes and (RT)-(q/d)PCR
                                assays from a multiple DNA sequence alignment. The workflow is
                                developed primarily for sequence variable RNA viruses, but it should also
                                be useful for other targets with high sequence variability"
                            ),
                            br(),
                            h5(tags$b("Instructions for use")),
                            h5("Add link to package vignette and manual"),
                            br(),
                            h5(tags$b("Citation")),
                            h5("S Persson et al., 2021, manuscript in preparation"),
                            br(),
                            h5(tags$b("R package and source code")),
                            uiOutput("codelink"),
                            br(),
                            h5(tags$b("Contact")),
                            h5("sofia.persson@slv.se")

                    ))
            ),

            # Target alignment =================================================

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

            # Consensus ========================================================

            tabItem(tabName = "consensus",
                    fluidRow(
                        column(width = 3, uiOutput("consensusSettings")),
                        column(width =  9, uiOutput("consensusOutput"))
                    )
            ),

            # Oligos ===========================================================

            tabItem(tabName = "oligos",
                    fluidRow(
                        column(width = 3, uiOutput("oligoSettings")),
                        column(width = 9,
                               uiOutput("oligoFilter"),
                               uiOutput("oligoOutput"),
                               uiOutput("oligoSelection")
                        ),
                    )
            ),

            # Assays ===========================================================

            tabItem(tabName = "assays",
                    fluidRow(
                        column(width = 3, uiOutput("assaySettings")),
                        column(width = 9,
                               uiOutput("assayFilter"),
                               uiOutput("assayOutput"),
                               uiOutput("assaySelection")
                        )
                    )
            )
        )
    )
)

# Server =======================================================================

server <- function(input, output) {

    # Data =====================================================================

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

    oligoCandidates <- eventReactive(input$getOligos, {
        req(consensus())
        oligos(consensus(),
               maxGapFrequency = input$maxGapFrequency,
               lengthPrimer = input$lengthPrimer,
               maxDegeneracyPrimer = input$maxDegeneracyPrimer,
               avoidThreeEndRunsPrimer = input$avoidThreeEndRunsPrimer,
               gcClampPrimer = input$gcClampPrimer,
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
        selectOligos(oligoCandidates(),
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

    selectedOligo <- reactive({
        req(oligoSelection())
        if (!is.null(input$table2_rows_selected)) {
            oligoSelection()[input$table2_rows_selected, ]
        } else {
            NULL
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
        req(oligoSelection())
        assays(oligoSelection(),
               length = input$length,
               tmDifferencePrimers = as.numeric(input$tmDifferencePrimers))
    })

    assaySelection <- reactive({
        req(assayCandidates())
        x <- assayCandidates()
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

    selectedAssay <- reactive({
        req(assaySelection())
        if (!is.null(input$table5_rows_selected)) {
            assaySelection()[input$table5_rows_selected, ]
        } else {
            NULL
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

    selectedAssayMatchList <- reactive({
        req(selectedAssayMatch())
        splitAssayToList(selectedAssayMatch())
    })

    # Links ====================================================================

    url <- a("GitHub", href = "https://github.com/sofpn/rprimer")
    output$codelink <- renderUI({
        tagList(url)
    })

    # UI output ================================================================


    #### Consensus settings ####

    output$consensusSettings <- renderUI({
        req(aln())
        box(width = 12, title = "Settings",
            numericInput(
                "ambiguityThreshold",
                h5("Threshold for an ambiguous base"),
                value = 0.05, min = 0, max = 0.2,
            ),
            h5("Region of interest"),
            numericInput(
                "roiFrom", h5("From"),
                min = 1, max = ncol(aln()) - 1, value = 1,

            ),
            numericInput(
                "roiTo", h5("To"),
                min = 2, max = ncol(aln()), value = ncol(aln())
            ),
            br(),
            actionButton("getConsensusProfile", "Get consensus profile")
        )
    })


    #### Consensus output ####

    output$consensusOutput <- renderUI({
        req(consensus())
        box(width = 12, title = "Consensus profile",
            tabBox(width = 12, title = "",
                   tabPanel(title = "Plot",
                            br(),
                            shinycssloaders::withSpinner(
                                plotOutput("plot1",
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
            )
        )
    })

    #### Oligo settings ####

    output$oligoSettings <- renderUI({
        req(consensus())
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
                choices = c(
                    "Ambiguous" = "ambiguous", "Mixed" = "mixed"
                    ),
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
            actionButton("getOligos", "Get oligos")
        )
    })

    #### Oligo filter #####

    output$oligoFilter <- renderUI({
        req(oligoCandidates())
        box(width = 12, title = "Filter",
            column(width = 12,
                   column(width = 2,
                          numericInput("fwdRegionFrom", h5("Fwd, from"),
                                       min = oligoCandidates()$roiStart[[1]],
                                       max = oligoCandidates()$roiEnd[[1]],
                                       value = oligoCandidates()$roiStart[[1]]
                          )
                   ),
                   column(width = 2,
                          numericInput("fwdRegionTo", h5("Fwd, to"),
                                       min = oligoCandidates()$roiStart[[1]],
                                       max = oligoCandidates()$roiEnd[[1]],
                                       value = oligoCandidates()$roiEnd[[1]]
                          )
                   ),
                   column(width = 2,
                          numericInput("revRegionFrom", h5("Rev, from"),
                                       min = oligoCandidates()$roiStart[[1]],
                                       max = oligoCandidates()$roiEnd[[1]],
                                       value = oligoCandidates()$roiStart[[1]]
                          )
                   ),
                   column(width = 2,
                          numericInput("revRegionTo", h5("Rev, to"),
                                       min = oligoCandidates()$roiStart[[1]],
                                       max = oligoCandidates()$roiEnd[[1]],
                                       value = oligoCandidates()$roiEnd[[1]]
                          )
                   ),
                   column(width = 2, uiOutput("prRegionFrom")),
                   column(width = 2, uiOutput("prRegionTo"))
            ),
            column(width = 12,
                   column(width = 4,
                          sliderInput("minOligoIdentityFwd",
                                      round = -4, step = 0.0001, ticks = FALSE,
                                      h5("Fwd, minimum identity"),
                                      min = round(min(
                                          oligoCandidates()$identity[
                                              oligoCandidates()$fwd &
                                                  oligoCandidates()$type == "primer"
                                          ],
                                          na.rm = TRUE), 4
                                      ),
                                      max = round(max(
                                          oligoCandidates()$identity[
                                              oligoCandidates()$fwd &
                                                  oligoCandidates()$type == "primer"
                                          ],
                                          na.rm = TRUE), 4
                                      ),
                                      value = round(min(
                                          oligoCandidates()$identity[
                                              oligoCandidates()$fwd &
                                                  oligoCandidates()$type == "primer"
                                          ],
                                          na.rm = TRUE), 4
                                      ),
                                      width = 200),
                          sliderInput("minOligoCoverageFwd",
                                      round = -4, step = 0.0001, ticks = FALSE,
                                      h5("Fwd, minimum coverage"),
                                      min = round(min(
                                          oligoCandidates()$coverage[
                                              oligoCandidates()$fwd &
                                                  oligoCandidates()$type == "primer"
                                              ],
                                          na.rm = TRUE), 4
                                      ),
                                      max = round(max(
                                          oligoCandidates()$coverage[
                                              oligoCandidates()$fwd &
                                                  oligoCandidates()$type == "primer"
                                          ],
                                          na.rm = TRUE), 4
                                      ),
                                      value = round(min(
                                          oligoCandidates()$coverage[
                                              oligoCandidates()$fwd &
                                                  oligoCandidates()$type == "primer"
                                          ], na.rm = TRUE), 4
                                      ),
                                      width = 200)
                   ),
                   column(width = 4,
                          sliderInput("minOligoIdentityRev",
                                      round = -4, step = 0.0001, ticks = FALSE,
                                      h5("Rev, minimum identity"),
                                      min = round(min(
                                          oligoCandidates()$identity[
                                              oligoCandidates()$rev &
                                                  oligoCandidates()$type == "primer"
                                              ],
                                          na.rm = TRUE), 4
                                      ),
                                      max = round(max(
                                          oligoCandidates()$identity[
                                              oligoCandidates()$rev &
                                                  oligoCandidates()$type == "primer"
                                          ],
                                          na.rm = TRUE), 4
                                      ),
                                      value = round(min(
                                          oligoCandidates()$identity[
                                              oligoCandidates()$rev &
                                                  oligoCandidates()$type == "primer"
                                          ],
                                          na.rm = TRUE), 4
                                      ),
                                      width = 200),
                          sliderInput("minOligoCoverageRev",
                                      round = -4, step = 0.0001, ticks = FALSE,
                                      h5("Rev, minimum coverage"),
                                      min = round(min(
                                          oligoCandidates()$coverage[
                                              oligoCandidates()$rev &
                                                  oligoCandidates()$type == "primer"
                                          ],
                                          na.rm = TRUE), 4
                                      ),
                                      max = round(max(
                                          oligoCandidates()$coverage[
                                              oligoCandidates()$rev &
                                                  oligoCandidates()$type == "primer"
                                          ],
                                          na.rm = TRUE), 4
                                      ),
                                      value = round(min(
                                          oligoCandidates()$coverage[
                                              oligoCandidates()$rev &
                                                  oligoCandidates()$type == "primer"
                                          ], na.rm = TRUE), 4
                                      ),
                                      width = 200)

                   ),
                   column(width = 4,
                          uiOutput("minOligoIdentityPr"),
                          uiOutput("minOligoCoveragePr")
                          )
            )
        )
    })

    output$prRegionFrom <- renderUI({
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

    output$minOligoIdentityPr <- renderUI({
        if (any(oligoCandidates()$type == "probe")) {
            sliderInput("minOligoIdentityPr",
                        round = -4, step = 0.0001, ticks = FALSE,
                        h5("Pr, minimum identity"),
                        min = round(min(
                            oligoCandidates()$identity[
                                oligoCandidates()$type == "probe"
                                ],
                            na.rm = TRUE), 4
                        ),
                        max = round(max(
                            oligoCandidates()$identity[
                                oligoCandidates()$type == "probe"
                            ],
                            na.rm = TRUE), 4
                        ),
                        value = round(min(
                            oligoCandidates()$identity[
                                oligoCandidates()$type == "probe"
                            ],
                            na.rm = TRUE), 4
                        ), width = 200)
        } else {
            NULL
        }

    })

    output$minOligoCoveragePr <- renderUI({
        if (any(oligoCandidates()$type == "probe")) {
            sliderInput("minOligoCoveragePr",
                        round = -4, step = 0.0001, ticks = FALSE,
                        h5("Pr, minimum coverage"),
                        min = round(min(
                            oligoCandidates()$coverage[
                                oligoCandidates()$type == "probe"
                            ],
                            na.rm = TRUE), 4),
                        max = round(max(
                            oligoCandidates()$coverage[
                                oligoCandidates()$type == "probe"
                            ],
                            na.rm = TRUE), 4),
                        value = round(min(
                            oligoCandidates()$coverage[
                                oligoCandidates()$type == "probe"
                            ],
                            na.rm = TRUE), 4),
                        width = 200)

        } else {
            NULL
        }

    })

    output$prRegionTo <- renderUI({
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

    #### Oligo output ####

    output$oligoOutput <- renderUI({
        req(oligoCandidates())
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
                            h5("Select a row for more details"),
                            br(),
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
                   )
            )
        )
    })

    #### Oligo selection ####

    output$oligoSelection <- renderUI({
        req(selectedOligo())
        box(width = 12, title = "Selection",
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
    })

    # Assays ===================================================================

    #### Assay settings ####

    output$assaySettings <- renderUI({
        req(oligoSelection())
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
            actionButton("getAssays", "Get assays")
        )
    })

    #### Assay filter ####

    output$assayFilter <- renderUI({
        req(assayCandidates())
        box(width = 12, title = "Filter",
            column(width = 12,
                   column(width = 2,
                          numericInput("assayRegionFrom", h5("From"),
                                       min = assayCandidates()$roiStart[[1]],
                                       max = assayCandidates()$roiEnd[[1]],
                                       value = assayCandidates()$roiStart[[1]]
                          )
                   ),
                   column(width = 2,
                          numericInput("assayRegionTo", h5("To"),
                                       min = assayCandidates()$roiStart[[1]],
                                       max = assayCandidates()$roiEnd[[1]],
                                       value = assayCandidates()$roiEnd[[1]]
                          )
                   )
            ),
            column(width = 12,
                   column(width = 3,
                          sliderInput("maxAssayScore",
                                      h5("Maximum score (lower is better)"),
                                      min = min(assayCandidates()$score, na.rm = TRUE),
                                      max = max(assayCandidates()$score, na.rm = TRUE),
                                      value = max(assayCandidates()$score, na.rm = TRUE),
                                      width = 200)
                   )

            )
        )
    })

    #### Assay output ####

    output$assayOutput <- renderUI({
        req(assayCandidates())
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
                            h5("Select a row for more details"),
                            br(),
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
                   )))
    })

    #### Assay selection ####

    output$assaySelection <- renderUI({
        req(selectedAssay())
        box(width = 12, title = "Selection",
            box(width = 12, title = "Overview",
                solidHeader = TRUE,
                downloadLink(
                    "download8", "Download summary table"
                ),
                br(),
                downloadLink(
                    "download9", "Download assay as fasta"
                ),
                br(),
                br(),
                DT::dataTableOutput("table6"),
                br(),
                br(),
                h5("Position in target alignment"),
                shinycssloaders::withSpinner(plotOutput(
                    "plotNtAssay",
                    height = plotHeight / 2,
                    width = "100%"
                ), color = "grey")
            ),
            box(width = 12, title = "Forward",
                solidHeader = TRUE,
                DT::dataTableOutput("tableFwd"),
                br(),
                h5("All sequence variants"),
                DT::dataTableOutput("table7"),
                br(),
                h5("Match information"),
                DT::dataTableOutput("tableFwdMatch"),
                br(),
                h5("Nucleotide distribution in target alignment"),
                shinycssloaders::withSpinner(plotOutput(
                    "plot7",
                    height = plotHeight / 2,
                    width = "100%"
                ), color = "grey"),
            ),
            box(width = 12, title = "Reverse",
                solidHeader = TRUE,
                DT::dataTableOutput("tableRev"),
                br(),
                h5("All sequence variants"),
                DT::dataTableOutput("table8"),
                br(),
                h5("Match information"),
                DT::dataTableOutput("tableRevMatch"),
                br(),
                h5("Nucleotide distribution in target alignment"),
                shinycssloaders::withSpinner(plotOutput(
                    "plot8",
                    height = plotHeight / 2,
                    width = "100%"
                ), color = "grey"),
            ),
            assayProbeInfo,
            box(width = 12, title = "Match details",
                solidHeader = TRUE,
                shinycssloaders::withSpinner(plotOutput(
                    "plot10",
                    height = plotHeight / 2,
                    width = "100%"
                ), color = "grey"),
                br(),
                h5("Forward"),
                br(),
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
                h5("Reverse"),
                br(),
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
                probeSelectionMatch
            )






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

    output$html3 <- renderText({
        req(!is.null(selectedOligoMatch()))
        c(
            "<b>Perfect match</b><br>",
            selectedOligoMatch()$idPerfectMatch[[1]])
    })

    output$html4 <- renderText({
        req(!is.null(selectedOligoMatch()))
        c(
            "<b>One mismatch</b><br>",
            selectedOligoMatch()$idOneMismatch[[1]])
    })

    output$html5 <- renderText({
        req(!is.null(selectedOligoMatch()))
        c(
            "<b>Two mismatches</b><br>",
            selectedOligoMatch()$idTwoMismatches[[1]])
    })

    output$html6 <- renderText({
        req(!is.null(selectedOligoMatch()))
        c(
            "<b>Three mismatches</b><br>",
            selectedOligoMatch()$idThreeMismatches[[1]])
    })

    output$html7 <- renderText({
        req(!is.null(selectedOligoMatch()))
        c(
            "<b>Four or more mismatches</b><br>",
            selectedOligoMatch()$idFourOrMoreMismatches[[1]])
    })

    output$html9 <- renderText({
        req(!is.null(selectedAssayMatch()))
        c(
            "<b>Perfect match</b><br>",
            selectedAssayMatch()$idPerfectMatchFwd[[1]])
    })

    output$html10 <- renderText({
        req(!is.null(selectedAssayMatch()))
        c(
            "<b>One mismatch</b><br>",
            selectedAssayMatch()$idOneMismatchFwd[[1]])
    })

    output$html11 <- renderText({
        req(!is.null(selectedAssayMatch()))
        c(
            "<b>Two mismatches</b><br>",
            selectedAssayMatch()$idTwoMismatchesFwd[[1]])
    })

    output$html12 <- renderText({
        req(!is.null(selectedAssayMatch()))
        c(
            "<b>Three mismatches</b><br>",
            selectedAssayMatch()$idThreeMismatchesFwd[[1]])
    })

    output$html13 <- renderText({
        req(!is.null(selectedAssayMatch()))
        c(
            "<b>Four or more mismatches</b><br>",
            selectedAssayMatch()$idFourOrMoreMismatchesFwd[[1]])
    })

    output$html14 <- renderText({
        req(!is.null(selectedAssayMatch()))
        c(
            "<b>Perfect match</b><br>",
            selectedAssayMatch()$idPerfectMatchRev[[1]])
    })

    output$html15 <- renderText({
        req(!is.null(selectedAssayMatch()))
        c(
            "<b>One mismatch</b><br>",
            selectedAssayMatch()$idOneMismatchRev[[1]])
    })

    output$html16 <- renderText({
        req(!is.null(selectedAssayMatch()))
        c(
            "<b>Two mismatches</b><br>",
            selectedAssayMatch()$idTwoMismatchesRev[[1]])
    })

    output$html17 <- renderText({
        req(!is.null(selectedAssayMatch()))
        c(
            "<b>Three mismatches</b><br>",
            selectedAssayMatch()$idThreeMismatchesRev[[1]])
    })

    output$html18 <- renderText({
        req(!is.null(selectedAssayMatch()))
        c(
            "<b>Four or more mismatches</b><br>",
            selectedAssayMatch()$idFourOrMoreMismatchesRev[[1]])
    })

    output$html19 <- renderText({
        req(any(grepl("Pr", names(selectedAssayMatch()))))
        c(
            "<b>Perfect match</b><br>",
            selectedAssayMatch()$idPerfectMatchPr[[1]])
    })

    output$html20 <- renderText({
        req(any(grepl("Pr", names(selectedAssayMatch()))))
        c(
            "<b>One mismatch</b><br>",
            selectedAssayMatch()$idOneMismatchPr[[1]])
    })

    output$html21 <- renderText({
        req(any(grepl("Pr", names(selectedAssayMatch()))))
        c(
            "<b>Two mismatches</b><br>",
            selectedAssayMatch()$idTwoMismatchesPr[[1]])
    })

    output$html22 <- renderText({
        req(any(grepl("Pr", names(selectedAssayMatch()))))
        c(
            "<b>Three mismatches</b><br>",
            selectedAssayMatch()$idThreeMismatchesPr[[1]])
    })

    output$html23 <- renderText({
        req(any(grepl("Pr", names(selectedAssayMatch()))))
        c(
            "<b>Four or more mismatches</b><br>",
            selectedAssayMatch()$idFourOrMoreMismatchesPr[[1]])
    })

    #### Plots ####

    output$plot1 <- renderPlot({
        plotData(consensus())
    })

    output$plot3 <- renderPlot({
        plotData(oligoSelection())
    })

    output$plot4 <- renderPlot({
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
        plotData(assaySelection())
    })

    output$plot7 <- renderPlot({
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
        if (is.na(selectedAssay()$length[[1]])) {
            NULL
        } else {
            from <- selectedAssay()$startPr
            to <- selectedAssay()$endPr
            rc <- ifelse(selectedAssay()$plusPr, FALSE, TRUE)
            plotData(consensus()[
                consensus()$position >= from & consensus()$position <= to,
            ], type = "nucleotide", rc = rc)
        }
    })

    output$plot10 <- renderPlot({
        req(!is.null(selectedAssayMatch()))
        plotData(selectedAssayMatch())
    })

    output$plotNtAssay <- renderPlot({
        if (is.na(selectedAssay()$length[[1]])) {
            NULL
        } else {
            from <- selectedAssay()$start
            to <- selectedAssay()$end
            plotData(consensus(), highlight = c(from, to))
        }
    })


    #### Tables ####

    output$table1 <- DT::renderDataTable({
        x <- roundDbls(as.data.frame(consensus()))
        names(x) <- c(
            "Position", "A", "C", "G", "T", "Other", "Gaps", "Majority",
            "Identity", "IUPAC", "Entropy", "Coverage"
        )
        x
    }, options = list(
        info = FALSE,
        searching = FALSE, paging = TRUE,
        scrollX = TRUE, autoWidth = FALSE,
        ordering = FALSE, scrollY = "300"
    ), rownames = FALSE, selection  = "none"
    )

    output$table2 <- DT::renderDataTable({
        if (is.na(oligoSelection()$length[[1]])) {
            NULL
        } else {
            x <- roundDbls(removeListColumns(as.data.frame(oligoSelection())))
            names(x) <- c(
                "Type", "Forward", "Reverse", "Start", "End", "Length",
                "IUPAC sequence", "IUPAC sequence, RC", "Identity",
                "Coverage", "Degeneracy", "GC content, mean",
                "GC content, range", "Tm, mean", "Tm, range", "Delta G, mean",
                "Delta G, range", "Design method", "Score", "ROI, start",
                "ROI, end"
            )
            x        }
    }, options = list(
        info = FALSE,
        searching = FALSE, paging = TRUE,
        scrollX = TRUE, autoWidth = TRUE,
        ordering = TRUE, scrollY = "300"
    ), rownames = FALSE,
    selection = list(mode = "single")
    )

    output$table3 <- DT::renderDataTable({
        if (is.na(selectedOligo()$length[[1]])) {
            NULL
        } else {
            x <- roundDbls(removeListColumns(as.data.frame(selectedOligo())))
            names(x) <- c(
                "Type", "Forward", "Reverse", "Start", "End", "Length",
                "IUPAC sequence", "IUPAC sequence, RC", "Identity",
                "Coverage", "Degeneracy", "GC content, mean",
                "GC content, range", "Tm, mean", "Tm, range", "Delta G, mean",
                "Delta G, range", "Design method", "Score", "ROI, start",
                "ROI, end"
            )
            x
        }
    }, options = list(
        info = FALSE,
        searching = FALSE, paging = FALSE,
        scrollX = TRUE, autoWidth = TRUE,
        ordering = FALSE
    ), rownames = FALSE, selection  = "none"
    )

    output$table4 <- DT::renderDataTable({
        if (is.na(selectedOligo()$length[[1]])) {
            NULL
        } else {
            x <- roundDbls(makeListTable(as.data.frame(selectedOligo())))
            names(x) <- c(
                "Sequence", "Sequence, RC", "GC content", "Tm", "Delta G"
            )
            x
        }
    }, options = list(
        info = FALSE,
        searching = FALSE, paging = FALSE,
        scrollX = TRUE, autoWidth = FALSE,
        ordering = FALSE
    ), rownames = FALSE, selection  = "none"
    )

    output$table5 <- DT::renderDataTable({
        if (is.na(assaySelection()$length[[1]])) {
            NULL
        } else {
            x <- assaySelection()
            x <- roundDbls(removeListColumns(as.data.frame(x)))
            if (any(grepl("Pr", names(x)))) {
                x$iupacSequencePr <- ifelse(x$plusPr, x$iupacSequencePr, NA)
                x$iupacSequenceRcPr <- ifelse(x$minusPr, x$iupacSequenceRcPr, NA)
            }
            names(x) <- if (any(grepl("Pr", names(x)))) {
                c(
                    "Start", "End", "Length", "Total degeneracy", "Score",
                    "Start, forward", "End, forward", "Length, forward",
                    "IUPAC sequence, forward", "Identity, forward",
                    "Coverage, forward", "Degeneracy, forward",
                    "GC content, mean, forward", "GC content, range, forward",
                    "Tm, mean, forward", "Tm, range, forward",
                    "Delta G, mean, forward", "Delta G, range, forward",
                    "Design method, forward",
                    "Start, reverse", "End, reverse", "Length, reverse",
                    "IUPAC sequence, reverse", "Identity, reverse",
                    "Coverage, reverse", "Degeneracy, reverse",
                    "GC content, mean, reverse", "GC content, range, reverse",
                    "Tm, mean, reverse", "Tm, range, reverse",
                    "Delta G, mean, reverse", "Delta G, range, reverse",
                    "Design method, reverse",
                    "Plus sense, probe", "Minus sense, probe",
                    "Start, probe", "End, probe", "Length, probe",
                    "IUPAC sequence, probe", "IUPAC sequence RC, probe",
                    "Identity, probe",
                    "Coverage, probe", "Degeneracy, probe",
                    "GC content, mean, probe", "GC content, range, probe",
                    "Tm, mean, probe", "Tm, range, probe",
                    "Delta G, mean, probe", "Delta G, range, probe",
                    "Design method, probe", "ROI, start", "ROI, end"
                )
            } else {
                c(
                    "Start", "End", "Length", "Total degeneracy", "Score",
                    "Start, forward", "End, forward", "Length, forward",
                    "IUPAC sequence, forward", "Identity, forward",
                    "Coverage, forward", "Degeneracy, forward",
                    "GC content, mean, forward", "GC content, range, forward",
                    "Tm, mean, forward", "Tm, range, forward",
                    "Delta G, mean, forward", "Delta G, range, forward",
                    "Design method, forward",
                    "Start, reverse", "End, reverse", "Length, reverse",
                    "IUPAC sequence, reverse", "Identity, reverse",
                    "Coverage, reverse", "Degeneracy, reverse",
                    "GC content, mean, reverse", "GC content, range, reverse",
                    "Tm, mean, reverse", "Tm, range, reverse",
                    "Delta G, mean, reverse", "Delta G, range, reverse",
                    "Design method, reverse",
                    "ROI, start", "ROI, end"
                )
            }
            x
        }
    }, options = list(
        info = FALSE,
        searching = FALSE, paging = TRUE,
        scrollX = TRUE, autoWidth = TRUE,
        ordering = TRUE, scrollY = "300"
    ), rownames = FALSE,
    selection = list(mode = "single")
    )

    output$table6 <- DT::renderDataTable({
        if (is.na(selectedAssay()$length[[1]])) {
            NULL
        } else {
            assayOverviewTable(selectedAssay())
        }
    }, options = list(
        info = FALSE,
        searching = FALSE, paging = FALSE,
        scrollX = TRUE, autoWidth = FALSE,
        ordering = FALSE
    ), rownames = FALSE, selection  = "none"
    )

    output$tableFwd <- DT::renderDataTable({
        x <- roundDbls(removeListColumns(selectedAssayList()[[1]]))
        names(x) <- c(
            "Start", "End", "Length", "IUPAC sequence", "Identity",
            "Coverage", "Degeneracy", "GC, mean", "GC, range", "Tm, mean",
            "Tm, range", "Delta G, mean", "Delta G, range", "Method"
            )
        x
    }, options = list(
        info = FALSE,
        searching = FALSE, paging = FALSE,
        scrollX = TRUE, autoWidth = FALSE,
        ordering = FALSE
    ), rownames = FALSE, selection  = "none"
    )

    output$tableFwdMatch <- DT::renderDataTable({
        x <- roundDbls(removeListColumns(selectedAssayMatchList()[[1]]))
        x <- x[ , -1]
        names(x) <- c(
            "Perfect match", "1 mismatch", "2 mismatches", "3 mismatches",
            "4 or more mismatches"
            )
        x
    }, options = list(
        info = FALSE,
        searching = FALSE, paging = FALSE,
        scrollX = TRUE, autoWidth = FALSE,
        ordering = FALSE
    ), rownames = FALSE, selection  = "none"
    )

    output$table7 <- DT::renderDataTable({
        x <- roundDbls(makeListTable(as.data.frame(selectedAssayList()[[1]])))
        names(x) <- c(
            "Sequence", "GC content", "Tm", "Delta G"
        )
        x
    }, options = list(
        info = FALSE,
        searching = FALSE, paging = FALSE,
        scrollX = TRUE, autoWidth = FALSE,
        ordering = FALSE
    ), rownames = FALSE, selection  = "none"
    )

    output$tableRev <- DT::renderDataTable({
        x <- roundDbls(removeListColumns(selectedAssayList()[[2]]))
        names(x) <- c(
            "Start", "End", "Length", "IUPAC sequence", "Identity",
            "Coverage", "Degeneracy", "GC, mean", "GC, range", "Tm, mean",
            "Tm, range", "Delta G, mean", "Delta G, range", "Method"
        )
        x
    }, options = list(
        info = FALSE,
        searching = FALSE, paging = FALSE,
        scrollX = TRUE, autoWidth = FALSE,
        ordering = FALSE
    ), rownames = FALSE, selection  = "none"
    )

    output$tableRevMatch <- DT::renderDataTable({
        x <- roundDbls(removeListColumns(selectedAssayMatchList()[[2]]))
        x <- x[ , -1]
        names(x) <- c(
            "Perfect match", "1 mismatch", "2 mismatches", "3 mismatches",
            "4 or more mismatches"
        )
        x
    }, options = list(
        info = FALSE,
        searching = FALSE, paging = FALSE,
        scrollX = TRUE, autoWidth = FALSE,
        ordering = FALSE
    ), rownames = FALSE, selection  = "none"
    )

    output$table8 <- DT::renderDataTable({
        x <- roundDbls(makeListTable(as.data.frame(selectedAssayList()[[2]])))
        names(x) <- c(
            "Sequence", "GC content", "Tm", "Delta G"
        )
        x
    }, options = list(
        info = FALSE,
        searching = FALSE, paging = FALSE,
        scrollX = TRUE, autoWidth = FALSE,
        ordering = FALSE
    ), rownames = FALSE, selection  = "none"
    )

    output$tablePr <- DT::renderDataTable({
        req(length(selectedAssayList()) == 3)
        x <- roundDbls(removeListColumns(selectedAssayList()[[3]]))
        x$iupacSequence <- ifelse(x$plus, x$iupacSequence, NA)
        x$iupacSequenceRc <- ifelse(x$minus, x$iupacSequenceRc, NA)
        x <- x[, c(-1, -2)]
        names(x) <- c(
            "Start", "End", "Length", "IUPAC sequence, plus",
            "Iupac sequence, minus",
            "Identity",
            "Coverage", "Degeneracy", "GC, mean", "GC, range", "Tm, mean",
            "Tm, range", "Delta G, mean", "Delta G, range", "Method"
        )
        x
    }, options = list(
        info = FALSE,
        searching = FALSE, paging = FALSE,
        scrollX = TRUE, autoWidth = TRUE,
        ordering = FALSE
    ), rownames = FALSE, selection  = "none"
    )

    output$tablePrMatch <- DT::renderDataTable({
        req(length(selectedAssayList()) == 3)
        x <- roundDbls(removeListColumns(selectedAssayMatchList()[[3]]))
        x <- x[ , -1]
        names(x) <- c(
            "Perfect match", "1 mismatch", "2 mismatches", "3 mismatches",
            "4 or more mismatches"
        )
        x
    }, options = list(
        info = FALSE,
        searching = FALSE, paging = FALSE,
        scrollX = TRUE, autoWidth = FALSE,
        ordering = FALSE
    ), rownames = FALSE, selection  = "none"
    )

    output$table9 <- DT::renderDataTable({
        req(length(selectedAssayList()) == 3)
        x <- roundDbls(makeListTable(as.data.frame(selectedAssayList()[[3]])))
        names(x) <- c(
            "Sequence, plus", "Sequence, minus", "GC content", "Tm", "Delta G"
        )
     #   if (!selectedAssayList()$plus) {
     #       x <- x[c("Sequence, minus", "GC content", "Tm", "Delta G")]
     #   }
     #   if (!selectedAssayList()$minus) {
     #       x <- x[c("Sequence, plus", "GC content", "Tm", "Delta G")]
     #   }
        x    }, options = list(
        info = FALSE,
        searching = FALSE, paging = FALSE,
        scrollX = TRUE, autoWidth = FALSE,
        ordering = FALSE
    ), rownames = FALSE, selection  = "none"
    )



    #### Download links ####

    output$download1 <- downloadHandler(
        filename <- function() {
            paste0("consensus_profile-", Sys.Date(), ".txt")
        },
        content <- function(file) {
            write.table(
                as.data.frame(consensus()), file,
                quote = FALSE, sep = "\t",
                row.names = FALSE
            )
        })

    output$download2 <- downloadHandler(
        filename <- function() {
            paste0("oligos-", Sys.Date(), ".txt")
        },
        content <- function(file) {
            write.table(
                as.data.frame(oligoSelection()), file,
                quote = FALSE, sep = "\t",
                row.names = FALSE
            )
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
            paste0("oligo-selection-", Sys.Date(), ".txt")
        },
        content <- function(file) {
            write.table(
                as.data.frame(selectedOligo()), file,
                quote = FALSE, sep = "\t",
                row.names = FALSE
            )
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
            paste0("assays-", Sys.Date(), ".txt")
        },
        content <- function(file) {
            write.table(
                as.data.frame(assaySelection()), file,
                quote = FALSE, sep = "\t",
                row.names = FALSE
            )
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
            paste0("assay-selection-", Sys.Date(), ".txt")
        },
        content <- function(file) {
            write.table(
                as.data.frame(selectedAssay()), file,
                quote = FALSE, sep = "\t",
                row.names = FALSE
            )
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
