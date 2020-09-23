context("sequence properties")

test_that("majority_consensus returns error when it should", {
  expect_error(majority_consensus(unclass(example_rprimer_profile)))
})

test_that("majority_consensus works", {
  cons <- majority_consensus(example_rprimer_profile)
  expect_equal(length(cons), ncol(example_rprimer_profile))
  expect_true(is.character(cons))
  expect_false(any(grepl("^acgt-", cons)))
})

test_that("asIUPAC returns an error when it should", {
  expect_error(asIUPAC(c("C,G", "G,T")))
})

test_that("asIUPAC works", {
  expect_true(is.na(asIUPAC("GC")))
  expect_equal(asIUPAC("C,G"), "S")
  expect_equal(asIUPAC("C,G,G"), "S")
  expect_equal(asIUPAC("C , G"), "S")
  expect_equal(asIUPAC("-"), "-")
})

test_that("iupac_consensus returns an error when it should", {
  expect_error(
    iupac_consensus(example_rprimer_profile, threshold = 3)
  )
  expect_error(
    iupac_consensus(example_rprimer_profile, threshold = NA)
  )
  expect_error(
    iupac_consensus(unclass(example_rprimer_profile))
  )
})

test_that("iupac_consensus works", {
  cons <- iupac_consensus(example_rprimer_profile, threshold = 0.2)
  expect_true(is.character(cons))
  expect_equal(length(cons), ncol(example_rprimer_profile))
  cons2 <- example_rprimer_profile[1:2, ]
  cons2 <- new_rprimer_profile(cons2)
  expect_warning(iupac_consensus(cons2, threshold = 0.2))
})

test_that("gap_frequency returns an error when it should", {
  expect_error(gap_frequency(unclass(example_rprimer_profile)))
})

test_that("gap_frequency works", {
  gaps <- gap_frequency(example_rprimer_profile)
  expect_true(is.double(gaps))
  expect_false(any(gaps < 0))
  expect_false(any(gaps > 1))
})

test_that("nucleotide_identity returns an error when it should", {
  expect_error(nucleotide_identity(unclass(example_rprimer_profile)))
})

test_that("nucleotide_identity works", {
  identity <- nucleotide_identity(example_rprimer_profile)
  expect_true(is.double(identity))
  expect_false(any(identity < 0))
  expect_false(any(identity > 1))
  expect_false(any(is.na(identity)))
  bases <- c("a", "c", "g", "t")
  s <- example_rprimer_profile[
    which(rownames(example_rprimer_profile) %in% bases),
    ]
  s <- new_rprimer_profile(s)
  expect_equal(nucleotide_identity(s), identity)
})

test_that("shannon_entropy returns an error when it should", {
  expect_error(shannon_entropy(unclass(example_rprimer_profile)))
})

test_that("shannon_entropy works", {
  entropy <- shannon_entropy(example_rprimer_profile)
  expect_true(is.double(entropy))
  expect_false(any(entropy < 0))
  expect_false(any(is.na(entropy)))
  bases <- c("a", "c", "g", "t")
  s <- example_rprimer_profile[
    which(rownames(example_rprimer_profile) %in% bases),
    ]
  s <- new_rprimer_profile(s)
  expect_equal(shannon_entropy(s), entropy)
})

test_that("sequence_properties returns an error when it should", {
  expect_error(
    sequence_properties(unclass(example_rprimer_profile))
  )
  expect_error(
    sequence_properties(example_rprimer_profile, iupac_threshold = 3)
  )
  expect_error(
    sequence_properties(example_rprimer_profile, iupac_threshold = NA)
  )
})

test_that("sequence_properties works", {
  properties <- sequence_properties(example_rprimer_profile)
  expect_s3_class(properties, class = "rprimer_properties")
  expect_equal(
    colnames(properties),
    c("position", "majority", "iupac", "gaps", "identity", "entropy")
  )
})
