welcomeUI <- function(id) {
    tagList(
        titlePanel("Welcome!"),
        h5(em("This application is under development")),
        h5(tags$b("Introduction")),
        h5(
            "This application provides tools for visualizing sequence
                                conservation and designing degenerate primers, probes and (RT)-(q/d)PCR
                                assays from a multiple DNA sequence alignment. The workflow is
                                developed primarily for sequence variable RNA viruses, but should also
                                be useful for other targets with high sequence variability"
        ),
        h5(tags$b("Instructions for use")),
        h5("Add link to package vignette and manual"),
        h5(tags$b("Citation")),
        h5("S Persson et al., 2021, manuscript in preparation"),
        h5(tags$b("R package and source code")),
        a("GitHub", href = "https://github.com/sofpn/rprimer"),
        h5(tags$b("Contact")),
        h5("sofia.persson@slv.se")
    )
}
