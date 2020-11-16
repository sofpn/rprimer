## Code to generate example data
library(magrittr)

exampleRprimerAlignment <- system.file(
  "extdata", "example_alignment.txt", package = "rprimer"
  ) %>%
    Biostrings::readDNAMultipleAlignment(., format = "fasta") %>%
    Biostrings::maskGaps(., min.fraction = 0.5, min.block.width = 1)

exampleRprimerProfile <- getConsensusProfile(
    exampleRprimerAlignment,
    iupacThreshold = 0.05
)

exampleRprimerOligo <- getOligos(exampleRprimerProfile,
                                 lengthPrimer = 18:22,
                                 maxGapFrequencyPrimer = 0.05,
                                 maxDegeneracyPrimer = 2,
                                 minEndIdentityPrimer = 1,
                                 gcRangePrimer = c(0.45, 0.65),
                                 tmRangePrimer = c(55, 65),
                                 concPrimer = 500,
                                 lengthProbe = 16:24,
                                 maxGapFrequencyProbe = 0.05,
                                 maxDegeneracyProbe = 4,
                                 avoid5EndGProbe = TRUE,
                                 gcRangeProbe = c(0.45, 0.65),
                                 tmRangeProbe = c(55, 70)
                                )

exampleRprimerAssay <- getAssays(
  exampleRprimerOligo, maxTmDifferencePrimers = 2,
  tmDifferencePrimersProbe = c(0, 10)
  )

save(exampleRprimerAlignment, file = "exampleRprimerAlignment.RData")
save(exampleRprimerProfile, file = "exampleRprimerProfile.RData")
save(exampleRprimerOligo, file = "exampleRprimerOligo.RData")
save(exampleRprimerAssay, file = "exampleRprimerAssay.RData")

tools::resaveRdaFiles(".")
