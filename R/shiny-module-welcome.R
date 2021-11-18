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
            "The design workflow in this application consist of five steps:"
        ),
        shiny::h5(shiny::tags$ul(
            shiny::tags$li("Generation of a consensus profile"),
            shiny::tags$li("Oligo design"),
            shiny::tags$li("Optional filtering of oligos"),
            shiny::tags$li("Assay design"),
            shiny::tags$li("Optional filtering of assays")
        )),
        shiny::h5(
            "Please see the package vingette or the package manual for
            further information on the different steps of the design procedure."
        ),
        shiny::hr()
    )
}
