welcomeUI <- function(id) {
    shiny::tagList(
        shiny::h5("Welcome!"),
        shiny::hr(),
        shiny::h6(
            "This application provides tools for visualizing sequence
                                conservation and designing degenerate primers, probes and (RT)-(q/d)PCR
                                assays from a multiple DNA sequence alignment."
        ),
        shiny::h6(
            "The workflow consists of five steps:"
        ),
        shiny::h6(shiny::tags$ul(
            shiny::tags$li("Generation of a consensus profile"),
            shiny::tags$li("Oligo design"),
            shiny::tags$li("Optional filtering of oligos"),
            shiny::tags$li("Assay design"),
            shiny::tags$li("Optional filtering of assays")
        )),
        shiny::h6(
            "Please see the package vingette or the package manual for
            further information on the different steps of the design procedure."
        ),
        shiny::hr()
    )
}
