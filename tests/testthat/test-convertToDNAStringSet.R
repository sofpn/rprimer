test_that("convertToDNAStringSetWorks", {

    data("exampleRprimerProfile")
    data("exampleRprimerOligo")
    data("exampleRprimerAssay")

    assay <- exampleRprimerAssay[1:3, ]
    test <- convertToDNAStringSet(assay)
    expect_s4_class(test, "DNAStringSet")
    test <- convertToDNAStringSet(assay, asRc = FALSE)
    expect_s4_class(test, "DNAStringSet")
    test <- convertToDNAStringSet(assay, asRc = TRUE, asIUPAC = FALSE)
    expect_true(all(grepl("_variant_", names(test))))

    noProbes <- oligos(
        exampleRprimerProfile, maxDegeneracyPrimer = 2, probe = FALSE
    )
##    noProbes <- assays(noProbes) ######################################################################


    oligo <- exampleRprimerOligo[1:10, ]
    test <- convertToDNAStringSet(oligo, asRc = TRUE)
    expect_s4_class(test, "DNAStringSet")
    expect_equal(length(test), 10)
    expect_true(all(grepl("_rc", names(test))))
    test <- convertToDNAStringSet(oligo, asRc = FALSE, asIUPAC = FALSE)
    expect_true(all(grepl("_variant_", names(test))))
    test <- convertToDNAStringSet(oligo, asRc = TRUE, asIUPAC = FALSE)
    expect_true(all(grepl("_rc", names(test))))
})
