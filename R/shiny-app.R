rprimerApp <- function() {

    ui <- shiny::fluidPage(

        shiny::tags$head(
            shiny::tags$style(
                shiny::HTML(".shiny-notification {
             position:fixed;
             top: calc(50%);
             left: calc(50%);
             }
             ")
            )
        ),
        shinyFeedback::useShinyFeedback(),
        theme = bslib::bs_theme(bootswatch = "yeti"),
        shiny::titlePanel(
            title = "Design degenerate oligos from a multiple DNA sequence alignment"
        ),
        shiny::tabsetPanel(
            id = "wizard",
            type = "hidden",
            shiny::tabPanel(
                "page1", "",
                welcomeUI("welcome"),
                shiny::br(),
                shiny::fluidRow(
                    shiny::column(width = 6),
                    shiny::column(
                        width = 6, align = "right",
                        nextButton(2)
                    )
                ),
            ),
            shiny::tabPanel(
                "page2", "",
                dataUI("data"),
                shiny::br(),
                shiny::fluidRow(
                    shiny::column(
                        width = 6,
                        previousButton(1)
                    ),
                    shiny::column(
                        width = 6, align = "right",
                        shiny::uiOutput("page23")
                    )
                )
            ),
            shiny::tabPanel(
                "page3", "",
                consensusUI("consensus"),
                shiny::br(),
                shiny::fluidRow(
                    shiny::column(
                        width = 6,
                        previousButton(2)
                    ),
                    shiny::column(
                        width = 6, align = "right",
                        shiny::uiOutput("page34")
                    )
                )
            ),
            shiny::tabPanel(
                "page4", "",
                oligoUI("oligo"),
                shiny::br(),
                shiny::fluidRow(
                    shiny::column(
                        width = 6,
                        previousButton(3)
                    ),
                    shiny::column(
                        width = 6, align = "right",
                        shiny::uiOutput("page45")
                    )
                )
            ),
            shiny::tabPanel(
                "page5", "",
                oligoFilterUI("oligoFilter"),
                shiny::br(),
                shiny::fluidRow(
                    shiny::column(
                        width = 6,
                        previousButton(4)
                    ),
                    shiny::column(
                        width = 6, align = "right",
                        nextButton(6)
                    )
                )
            ),
            shiny::tabPanel(
                "page6", "",
                assayUI("assay"),
                shiny::br(),
                shiny::fluidRow(
                    shiny::column(
                        width = 6,
                        previousButton(5)
                    ),
                    shiny::column(
                        width = 6, align = "right",
                        shiny::uiOutput("page67")
                    )
                )
            ),
            shiny::tabPanel(
                "page7", "",
                assayFilterUI("assayFilter"),
                shiny::br(),
                shiny::fluidRow(
                    shiny::column(
                        width = 6,
                        previousButton(6)
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

        output$page23 <- shiny::renderUI({
            shiny::req(is(alignment$data(), "DNAMultipleAlignment"))
            nextButton(3)
        })

        output$page34 <- shiny::renderUI({
            shiny::req(is(consensus$data(), "RprimerProfile"))
            nextButton(4)
        })

        output$page45 <- shiny::renderUI({
            shiny::req(is(oligo$data(), "RprimerOligo"))
            nextButton(5)
        })

        output$page67 <- shiny::renderUI({
            shiny::req(is(assay$data(), "RprimerAssay"))
            nextButton(7)
        })

        shiny::observeEvent(input$page12, switchPage(2))
        shiny::observeEvent(input$page21, switchPage(1))
        shiny::observeEvent(input$page23, switchPage(3))
        shiny::observeEvent(input$page32, switchPage(2))
        shiny::observeEvent(input$page34, switchPage(4))
        shiny::observeEvent(input$page43, switchPage(3))
        shiny::observeEvent(input$page45, switchPage(5))
        shiny::observeEvent(input$page54, switchPage(4))
        shiny::observeEvent(input$page56, switchPage(6))
        shiny::observeEvent(input$page65, switchPage(5))
        shiny::observeEvent(input$page67, switchPage(7))
        shiny::observeEvent(input$page76, switchPage(6))
    }

    shiny::shinyApp(ui, server)

}

# Helpers ======================================================================

switchPage <- function(i) {
    shiny::updateTabsetPanel(inputId = "wizard", selected = paste0("page", i))
}

previousButton <- function(i) {
    id <- paste0("page", i + 1, i)
    shiny::actionButton(
        id, "Previous",
        icon = shiny::icon("chevron-left")
    )
}

nextButton <- function(i) {
    id <- paste0("page", i - 1, i)
    shiny::actionButton(
        id,
        label = shiny::div(
            "Next", shiny::icon("chevron-right"),
            class = "btn btn-primary"
        )
    )
}
