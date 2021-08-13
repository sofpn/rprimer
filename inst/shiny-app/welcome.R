welcomeUI <- function(id) {
    tagList(
        titlePanel("Welcome!"),
        h5(tags$b("Overview")),
        h5(
            "This application provides tools for visualizing sequence
                                conservation and designing degenerate primers, probes and (RT)-(q/d)PCR
                                assays from a multiple DNA sequence alignment. The workflow is
                                developed primarily for sequence variable RNA viruses, but it should also
                                be useful for other targets with high sequence variability"
        ),
        h5(tags$b("Citation")),
        h5("S Persson et al., 2021, manuscript in preparation"),
        tags$h5(tags$b(tags$a(href = "https://github.com/sofpn/rprimer", "View source", icon("github"))))
    )
}
