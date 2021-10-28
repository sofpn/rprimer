test_that("Data import and workflow works.", {
    infile <- system.file("extdata", "example_alignment.txt", package = "rprimer")
    aln <- Biostrings::readDNAMultipleAlignment(infile, format = "fasta")
    consensus <- consensusProfile(aln, 0.1)
    ols <- designOligos(consensus[1:2000, ])
    expect_s4_class(ols, "RprimerOligo")
})
