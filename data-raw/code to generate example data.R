## Code to generate example data

exampleRprimerAlignment <- system.file(
    "extdata", "example_alignment.txt",
    package = "rprimer"
)

exampleRprimerAlignment <- Biostrings::readDNAMultipleAlignment(
    exampleRprimerAlignment,
    format = "fasta"
)

exampleRprimerProfile <- consensusProfile(exampleRprimerAlignment, 0.05)

exampleRprimerOligo <- oligos(exampleRprimerProfile)

exampleRprimerAssay <- assays(exampleRprimerOligo)

exampleRprimerMatchOligo <- checkMatch(
    exampleRprimerOligo[1:10, ],
    target = exampleRprimerAlignment
)

exampleRprimerMatchAssay <- checkMatch(
    exampleRprimerAssay[1:5, ],
    target = exampleRprimerAlignment
)

save(exampleRprimerAlignment, file = "exampleRprimerAlignment.RData")

save(exampleRprimerProfile, file = "exampleRprimerProfile.RData")

save(exampleRprimerOligo, file = "exampleRprimerOligo.RData")

save(exampleRprimerAssay, file = "exampleRprimerAssay.RData")

save(exampleRprimerMatchOligo, file = "exampleRprimerMatchOligo.RData")

save(exampleRprimerMatchAssay, file = "exampleRprimerMatchAssay.RData")

tools::resaveRdaFiles(".")
