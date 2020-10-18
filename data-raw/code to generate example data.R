## Code to generate example data
exampleRprimerAlignment <- "exampleAlignment.txt" %>%
    Biostrings::readDNAMultipleAlignment(., format = "fasta") %>%
    Biostrings::maskGaps(., min.fraction = 0.5, min.block.width = 1)

exampleRprimerProfile <- getConsensusProfile(
  exampleRprimerAlignment, iupacThreshold = 0.05
)

save(exampleRprimerAlignment, file = "exampleRprimerAlignment.RData")
save(exampleRprimerProfile, file = "exampleRprimerProfile.RData")

tools::resaveRdaFiles(".")
