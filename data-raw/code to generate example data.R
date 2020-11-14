## Code to generate example data
library(magrittr)

exampleRprimerAlignment <- system.file("extdata", "example_alignment.txt", package = "rprimer") %>%
    Biostrings::readDNAMultipleAlignment(., format = "fasta") %>%
    Biostrings::maskGaps(., min.fraction = 0.5, min.block.width = 1)

exampleRprimerProfile <- getConsensusProfile(
    exampleRprimerAlignment,
    iupacThreshold = 0.05
)

exampleRprimerOligo <- getOligos(exampleRprimerProfile)

exampleRprimerAssay <- getAssays(exampleRprimerOligo)

save(exampleRprimerAlignment, file = "exampleRprimerAlignment.RData")
save(exampleRprimerProfile, file = "exampleRprimerProfile.RData")
save(exampleRprimerOligo, file = "exampleRprimerOligo.RData")
save(exampleRprimerAssay, file = "exampleRprimerAssay.RData")

tools::resaveRdaFiles(".")
