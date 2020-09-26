## Code to generate example data
exampleAlignment <- "exampleAlignment.txt" %>%
  Biostrings::readDNAMultipleAlignment(., format = "fasta") %>%
  Biostrings::maskGaps(., min.fraction = 0.5, min.block.width = 1)

exampleRprimerProfile <- getAlignmentProfile(exampleAlignment)
exampleRprimerProperties <- getAlignmentProperties(
  exampleRprimerProfile, iupacThreshold = 0.1
)

save(exampleRprimerProfile, file = "exampleRprimerProfile.RData")
save(exampleRprimerProperties, file = "exampleRprimerProperties.RData")

tools::resaveRdaFiles(".")
