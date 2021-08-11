library(shiny)
library(bslib)
library(rprimer)

source("utils.R")
source("welcome.R")
source("data.R")
source("consensus.R")
source("oligo.R")
source("oligo_filter.R")
source("assay.R")
source("assay_filter.R")

ui <- fluidPage(

    tags$head(
        tags$style(
            HTML(".shiny-notification {
             position:fixed;
             top: calc(50%);
             left: calc(50%);
             }
             ")
        )
    ),
    shinyFeedback::useShinyFeedback(),
    theme = bs_theme(bootswatch = "yeti"),
    titlePanel(
        title = "Design degenerate oligos from a multiple DNA sequence alignment"
    ),
    tabsetPanel(
        id = "wizard",
        type = "hidden",
        tabPanel(
            "page1", "",
            welcomeUI("welcome"),
            br(),
            fluidRow(
                column(width = 6),
                column(
                    width = 6, align = "right",
                    actionButton(
                        "page12",
                        label = div("Next", icon("chevron-right"))
                    )
                )
            ),
        ),
        tabPanel(
            "page2", "",
            dataUI("data"),
            br(),
            fluidRow(
                column(
                    width = 6,
                    actionButton(
                        "page21", "Previous",
                        icon = icon("chevron-left")
                    )
                ),
                column(
                    width = 6, align = "right",
                    uiOutput("page23")
                )
            )
        ),
        tabPanel(
            "page3", "",
            consensusUI("consensus"),
            br(),
            fluidRow(
                column(
                    width = 6,
                    actionButton(
                        "page32", "Previous",
                        icon = icon("chevron-left")
                    )
                ),
                column(
                    width = 6, align = "right",
                    uiOutput("page34")
                )
            )
        ),
        tabPanel(
            "page4", "",
            oligoUI("oligo"),
            br(),
            fluidRow(
                column(
                    width = 6,
                    actionButton(
                        "page43", "Previous",
                        icon = icon("chevron-left")
                    )
                ),
                column(
                    width = 6, align = "right",
                    uiOutput("page45")
                )
            )
        ),
        tabPanel(
            "page5", "",
            oligoFilterUI("oligoFilter"),
            br(),
            fluidRow(
                column(
                    width = 6,
                    actionButton(
                        "page54", "Previous",
                        icon = icon("chevron-left")
                    )
                ),
                column(
                    width = 6, align = "right",
                    uiOutput("page56")
                )
            )
        ),
        tabPanel(
            "page6", "",
            assayUI("assay"),
            br(),
            fluidRow(
                column(
                    width = 6,
                    actionButton(
                        "page65", "Previous",
                        icon = icon("chevron-left")
                    )
                ),
                column(
                    width = 6, align = "right",
                    uiOutput("page67")
                )
            )
        ),
        tabPanel(
            "page7", "",
            assayFilterUI("assayFilter"),
            br(),
            fluidRow(
                column(
                    width = 6,
                    actionButton(
                        "page76", "Previous",
                        icon = icon("chevron-left")
                    )
                )
            )
        )
    )
)

server <- function(input, output, session) {
    alignment <- dataServer("data")

    consensus <- consensusServer("consensus", alignment = alignment$data)

    oligo <- oligoServer(
        "oligo",
        alignment = alignment$data, consensus = consensus$data
    )

    oligoFilter <- oligoFilterServer(
        "oligoFilter",
        alignment = alignment$data, consensus = consensus$data,
        oligos = oligo$data
    )

    assay <- assayServer(
        "assay",
        alignment = alignment$data, consensus = consensus$data,
        oligos = oligoFilter$data
    )

    assayFilterServer(
        "assayFilter",
        alignment = alignment$data, consensus = consensus$data,
        assays = assay$data
    )

    output$page23 <- renderUI({
        req(is(alignment$data(), "DNAMultipleAlignment"))
        actionButton("page23", label = div("Next", icon("chevron-right")))
    })

    output$page34 <- renderUI({
        (is(consensus$data(), "RprimerProfile"))
        actionButton("page34", label = div("Next", icon("chevron-right")))
    })

    output$page45 <- renderUI({
        req(is(oligo$data(), "RprimerOligo"))
        actionButton("page45", label = div("Next", icon("chevron-right")))
    })

    output$page56 <- renderUI({
        req(is(oligoFilter$data(), "RprimerOligo"))
        actionButton("page56", label = div("Next", icon("chevron-right")))
    })

    output$page67 <- renderUI({
        req(is(assay$data(), "RprimerAssay"))
        actionButton("page67", label = div("Next", icon("chevron-right")))
    })

    switchPage <- function(i) {
        updateTabsetPanel(inputId = "wizard", selected = paste0("page", i))
    }

    observeEvent(input$page12, switchPage(2))
    observeEvent(input$page21, switchPage(1))
    observeEvent(input$page23, switchPage(3))
    observeEvent(input$page32, switchPage(2))
    observeEvent(input$page34, switchPage(4))
    observeEvent(input$page43, switchPage(3))
    observeEvent(input$page45, switchPage(5))
    observeEvent(input$page54, switchPage(4))
    observeEvent(input$page56, switchPage(6))
    observeEvent(input$page65, switchPage(5))
    observeEvent(input$page67, switchPage(7))
    observeEvent(input$page76, switchPage(6))
}

shinyApp(ui, server)
