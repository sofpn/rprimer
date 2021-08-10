options(shiny.maxRequestSize = 90 * 1024^2)

dataUI <- function(id) {
    ns <- NS(id)

    tagList(
        titlePanel("File import"),
        br(),
        sidebarLayout(
            sidebarPanel(
                radioButtons(
                    ns("dataSelection"),
                    label = NULL,
                    choices = c(
                        "Upload alignment", "Use example data"
                    ), selected = "Upload alignment"
                ),
                conditionalPanel(
                    ns = NS(id),
                    condition = "input.dataSelection == 'Upload alignment'",
                    radioButtons(
                        ns("filetype"), h5("Format"),
                        choices = c("Fasta" = "fasta", "Clustal" = "clustal"),
                        selected = "fasta"
                    ),
                    fileInput(
                        ns("file"), "",
                        multiple = FALSE,
                        accept = c("text", "fasta")
                    ),
                    htmlOutput(ns("fileName"))
                ),
                conditionalPanel(
                    ns = NS(id),
                    condition = "input.dataSelection == 'Use example data'",
                    h5("Hepatitis E virus alignment")
                )
            ),
            mainPanel(
                htmlOutput(ns("nSequences")),
                htmlOutput(ns("alnLength"))
            )
        )
    )
}

dataServer <- function(id) {
    moduleServer(id, function(input, output, session) {
        data("exampleRprimerAlignment")

        aln <- reactive({
            if (input$dataSelection == "Upload alignment") {
                req(input$file)
                tryCatch(
                    {
                        Biostrings::readDNAMultipleAlignment(
                            input$file$datapath,
                            format = input$filetype
                        )
                    },
                    error = function(cond) {
                        showNotification(
                            "Failed to upload file correctly.\n
                        Please check that the file format is correct.",
                            type = "error",
                            duration = NULL
                        )
                    },
                    silent = TRUE
                )
            } else if (input$dataSelection == "Use example data") {
                exampleRprimerAlignment
            }
        })

        output$fileName <- renderText({
            req(is(aln(), "DNAMultipleAlignment"))
            if (input$dataSelection == "Upload alignment") {
                paste0(
                    "'", input$file[[1]], "' was successfully uploaded! "
                )
            } else {
                NULL
            }
        })

        output$nSequences <- renderText({
            req(is(aln(), "DNAMultipleAlignment"))
            paste("Number of sequences:", nrow(aln()))
        })

        output$alnLength <- renderText({
            req(is(aln(), "DNAMultipleAlignment"))
            paste("Alignment length:", ncol(aln()))
        })

        list(data = reactive(aln()))
    })
}
