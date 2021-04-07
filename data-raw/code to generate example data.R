## Code to generate example data

exampleRprimerAlignment <- system.file(
    "extdata", "example_alignment.txt",
    package = "rprimer"
)

exampleRprimerAlignment <- Biostrings::readDNAMultipleAlignment(
    exampleRprimerAlignment,
    format = "fasta"
)

#exampleRprimerAlignment <- Biostrings::maskGaps(
#    exampleRprimerAlignment,
#    min.fraction = 0.5, min.block.width = 1
#)

exampleRprimerProfile <- consensusProfile(exampleRprimerAlignment, 0.05)
exampleRprimerOligo <- oligos(exampleRprimerProfile)
exampleRprimerAssay <- assays(exampleRprimerOligo)

save(exampleRprimerAlignment, file = "exampleRprimerAlignment.RData")
save(exampleRprimerProfile, file = "exampleRprimerProfile.RData")
save(exampleRprimerOligo, file = "exampleRprimerOligo.RData")
save(exampleRprimerAssay, file = "exampleRprimerAssay.RData")

tools::resaveRdaFiles(".")
