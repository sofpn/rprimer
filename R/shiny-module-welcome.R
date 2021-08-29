welcomeUI <- function(id) {
    shiny::tagList(
        shiny::h4("Welcome!"),
        shiny::hr(),
        shiny::h5(shiny::tags$b("Overview")),
        shiny::h5(
            "This application provides tools for visualizing sequence
                                conservation and designing degenerate primers, probes and (RT)-(q/d)PCR
                                assays from a multiple DNA sequence alignment."
        ),
        shiny::h5(
            "The design workflow consist of five steps:"
            ),
        shiny::h5(shiny::tags$ul(
            shiny::tags$li("Generation of a consensus profile"),
            shiny::tags$li("Oligo design"),
            shiny::tags$li("Optional filtering of oligos"),
            shiny::tags$li("Assay design"),
            shiny::tags$li("Optional filtering of assays")
        )),
        shiny::h5(shiny::tags$b("Citation")),
        shiny::h5("Persson et al., 2021, manuscript in preparation"),
        #tags$h5(tags$b(tags$a(href = "https://github.com/sofpn/rprimer", "View source", icon("github")))),
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
