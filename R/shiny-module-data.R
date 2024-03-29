dataUI <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h5("File import"),
        shiny::hr(),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::radioButtons(
                    ns("dataSelection"),
                    shiny::h6(shiny::tags$b("Option")),
                    choices = c(
                        "Upload alignment", "Use example data"
                    ), selected = "Upload alignment"
                ),
                shiny::conditionalPanel(
                    ns = shiny::NS(id),
                    condition = "input.dataSelection == 'Upload alignment'",
                    shiny::radioButtons(
                        ns("filetype"), shiny::h6(shiny::tags$b("File format")),
                        choices = c("Fasta" = "fasta", "Clustal" = "clustal"),
                        selected = "fasta"
                    ),
                    shiny::fileInput(
                        ns("file"), "",
                        multiple = FALSE,
                        accept = c("text", "fasta")
                    ),
                    shiny::htmlOutput(ns("fileName"))
                ),
                shiny::conditionalPanel(
                    ns = shiny::NS(id),
                    condition = "input.dataSelection == 'Use example data'",
                    shiny::h6("Hepatitis E virus alignment")
                )
            ),
            shiny::mainPanel(
                shiny::htmlOutput(ns("nSequences")),
                shiny::htmlOutput(ns("alnLength"))
            )
        ),
        shiny::hr()
    )
}

dataServer <- function(id) {
    shiny::moduleServer(id, function(input, output, session) {
        options(shiny.maxRequestSize=30*1024^2)

        aln <- shiny::reactive({
            if (input$dataSelection == "Upload alignment") {
                shiny::req(input$file)
                tryCatch(
                    {
                        readDNAMultipleAlignment(
                            input$file$datapath,
                            format = input$filetype
                        )
                    },
                    error = function(cond) {
                        shiny::showNotification(
                            "Failed to import sequence data from the file.\n
                        Please check that the file format is correct.",
                            type = "error",
                            duration = NULL
                        )
                    },
                    silent = TRUE
                )
            } else if (input$dataSelection == "Use example data") {
                utils::data("exampleRprimerAlignment", package = "rprimer")
                exampleRprimerAlignment
            }
        })

        output$fileName <- shiny::renderText({
            shiny::req(is(aln(), "DNAMultipleAlignment"))
            if (input$dataSelection == "Upload alignment") {
                paste0(
                    "'", input$file[[1]], "' was successfully uploaded! "
                )
            } else {
                NULL
            }
        })

        output$nSequences <- shiny::renderText({
            shiny::req(is(aln(), "DNAMultipleAlignment"))
            paste("Number of sequences:", nrow(aln()))
        })

        output$alnLength <- shiny::renderText({
            shiny::req(is(aln(), "DNAMultipleAlignment"))
            paste("Alignment width:", ncol(aln()))
        })

        list(data = shiny::reactive(aln()))
    })
}
