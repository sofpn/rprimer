# Import alignment and remove positions with high gap frequency
myAlignment <- "example_alignment.txt" %>%
  Biostrings::readDNAMultipleAlignment(., format = "fasta") %>%
  Biostrings::maskGaps(., min.fraction = 0.5, min.block.width = 1)

# Select region of interest
Biostrings::colmask(myAlignment, invert = TRUE) <- IRanges::IRanges(start = 1000, end = 5000)

#
myProfile <- Biostrings::consensusMatrix(myAlignment, as.prob = TRUE) %>%
  RprimerProfile(.)
