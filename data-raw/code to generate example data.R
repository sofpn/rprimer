## Code to generate example data

exampleRprimerAlignment <- system.file(
    "extdata", "example_alignment.txt",
    package = "rprimer"
)

exampleRprimerAlignment <- Biostrings::readDNAMultipleAlignment(
    exampleRprimerAlignment,
    format = "fasta"
)

exampleRprimerAlignment <- Biostrings::maskGaps(
    exampleRprimerAlignment,
    min.fraction = 0.5, min.block.width = 1
)

exampleRprimerProfile <- consensusProfile(
    exampleRprimerAlignment, 0.05
)

exampleRprimerOligo <- oligos(exampleRprimerProfile,
    lengthPrimer = 18:22,
    maxGapFrequency = 0.05,
    maxDegeneracyPrimer = 2,
    minEndIdentityPrimer = 0.9,
    gcRangePrimer = c(0.45, 0.65),
    tmRangePrimer = c(55, 65),
    concPrimer = 500,
    lengthProbe = 16:24,
    maxDegeneracyProbe = 4,
    avoidFiveEndGProbe = TRUE,
    gcRangeProbe = c(0.45, 0.65),
    tmRangeProbe = c(55, 70)
)

exampleRprimerAssay <- assays(
    exampleRprimerOligo,
    tmDiffPrimers = 2,
    tmDiffPrimersProbe = c(0, 10)
)

save(exampleRprimerAlignment, file = "exampleRprimerAlignment.RData")
save(exampleRprimerProfile, file = "exampleRprimerProfile.RData")
save(exampleRprimerOligo, file = "exampleRprimerOligo.RData")
save(exampleRprimerAssay, file = "exampleRprimerAssay.RData")

tools::resaveRdaFiles(".")
