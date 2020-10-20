## Code to generate example data
exampleRprimerAlignment <- "exampleAlignment.txt" %>%
    Biostrings::readDNAMultipleAlignment(., format = "fasta") %>%
    Biostrings::maskGaps(., min.fraction = 0.5, min.block.width = 1)

exampleRprimerProfile <- getConsensusProfile(
    exampleRprimerAlignment,
    iupacThreshold = 0.05
)

exampleRprimerOligo <- getOligos(exampleRprimerProfile)

exampleRprimerAssay <- getAssays(exampleRprimerOligo, exampleRprimerOligo)

save(exampleRprimerAlignment, file = "exampleRprimerAlignment.RData")
save(exampleRprimerProfile, file = "exampleRprimerProfile.RData")
save(exampleRprimerOligo, file = "exampleRprimerOligo.RData")
save(exampleRprimerAssay, file = "exampleRprimerAssay.RData")

tools::resaveRdaFiles(".")
