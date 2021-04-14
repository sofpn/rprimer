test_that("convertToDNAStringSetWorks", {
    data("exampleRprimerAssay")
    assay <- exampleRprimerAssay[1, ]
    assay <- convertToDNAStringSet(assay)
    expect_s4_class(assay, "DNAStringSet")

    assayRc <- exampleRprimerAssay[1, ]
    assayRc <- convertToDNAStringSet(assayRc, asRc = FALSE)
    expect_s4_class(assayRc, "DNAStringSet")

    data("exampleRprimerOligo")
    oligo <- exampleRprimerOligo[1:10, ]
    oligo <- convertToDNAStringSet(oligo, asRc = FALSE)
    expect_s4_class(oligo, "DNAStringSet")
})
