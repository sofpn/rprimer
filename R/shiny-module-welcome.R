welcomeUI <- function(id) {
    shiny::tagList(
        shiny::titlePanel("Welcome!"),
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
    )
}
