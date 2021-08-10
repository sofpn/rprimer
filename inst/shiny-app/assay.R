assayUI <- function(id) {
    ns <- NS(id)

    tagList(
        titlePanel("Step 4/5: Design assays"),
        br(),
        sidebarLayout(
            sidebarPanel(),
            mainPanel()
        )
    )
}

assayServer <- function(id, alignment, consensus, oligos) {
    moduleServer(id, function(input, output, session) {

    })
}
