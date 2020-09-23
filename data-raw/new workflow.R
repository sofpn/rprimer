# Import alignment and remove positions with high gap frequency
myAlignment <- "example_alignment.txt" %>%
  Biostrings::readDNAMultipleAlignment(., format = "fasta") %>%
  Biostrings::maskGaps(., min.fraction = 0.5, min.block.width = 1)

# Select region of interest
Biostrings::colmask(myAlignment, invert = TRUE) <- IRanges::IRanges(start = 1000, end = 5000)

#
myAlignmentProfile <- Biostrings::consensusMatrix(myAlignment, as.prob = TRUE) #%>%

x <- myAlignmentProfile

# Remove columns with NA
y <- x[, colSums(!is.na(x)) > 0]


# allClasses document
# vingette
# Rprimer
sessionInfo()
