assayFilterUI <- function(id) {
    ns <- NS(id)

    tagList(
        titlePanel("Step 5/5: Filter assays"),
        br(),
        sidebarLayout(
            sidebarPanel(),
            mainPanel()
        )
    )
}

assayFilterServer <- function(id, alignment, consensus, assays) {
    moduleServer(id, function(input, output, session) {

    })
}
