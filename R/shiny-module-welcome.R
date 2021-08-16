welcomeUI <- function(id) {
    shiny::tagList(
        shiny::h4("Welcome!"),
        shiny::hr(),
        shiny::h5(shiny::tags$b("Overview")),
        shiny::h5(
            "This application provides tools for visualizing sequence
                                conservation and designing degenerate primers, probes and (RT)-(q/d)PCR
                                assays from a multiple DNA sequence alignment. The workflow is
                                developed primarily for sequence variable RNA viruses, but it should also
                                be useful for other targets with high sequence variability"
        ),
        shiny::h5(shiny::tags$b("Citation")),
        shiny::h5("S Persson et al., 2021, manuscript in preparation"),
        shiny::hr()
    )
}

## Module app for testing ======================================================

# welcomeApp <- function() {
#    ui <- fluidPage(
#        welcomeUI("id")
#    )
#    server <- function(input, output, session) {}
#    shinyApp(ui, server)
# }
