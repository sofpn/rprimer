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
            oligo = oligo$data
        )

        assay <- assayServer(
            "assay",
            alignment = alignment$data, consensus = consensus$data,
            oligo = oligoFilter$data
        )

        assayFilterServer(
            "assayFilter",
            alignment = alignment$data, consensus = consensus$data,
            assay = assay$data
        )

        output$page23 <- shiny::renderUI({
            shiny::req(alignment$data())
            nextButton(3)
        })

        output$page34 <- shiny::renderUI({
            shiny::req(consensus$data())
            nextButton(4)
        })

        output$page45 <- shiny::renderUI({
            shiny::req(oligo$data())
            nextButton(5)
        })

        output$page56 <- shiny::renderUI({
            shiny::req(oligoFilter$data())
            nextButton(6)
        })

        output$page67 <- shiny::renderUI({
            shiny::req(assay$data())
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
