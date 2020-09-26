## Code to generate example data
exampleRprimerAlignment <- "exampleAlignment.txt" %>%
  Biostrings::readDNAMultipleAlignment(., format = "fasta") %>%
  Biostrings::maskGaps(., min.fraction = 0.5, min.block.width = 1)

exampleRprimerProfile <- getAlignmentProfile(exampleRprimerAlignment)
exampleRprimerProperties <- getAlignmentProperties(
  exampleRprimerProfile, iupacThreshold = 0.1
)

save(exampleRprimerAlignment, file = "exampleRprimerAlignment.RData")
save(exampleRprimerProfile, file = "exampleRprimerProfile.RData")
save(exampleRprimerProperties, file = "exampleRprimerProperties.RData")

tools::resaveRdaFiles(".")
