context("import alignment")

test_that("read_fasta_alignment stops when it should", {
  expect_error(read_fasta_alignment(invalid_filename))
  expect_error(read_fasta_alignment("file_not_found.txt"))
  expect_error(read_fasta_alignment("non_fasta_format.txt"))
  expect_error(read_fasta_alignment("not_unique_sequence_names.txt"))
  expect_error(read_fasta_alignment("different_sequence_lengths.txt"))
  expect_error(read_fasta_alignment("invalid_sequence_character.txt"))
  expect_error(read_fasta_alignment("too_short_alignment.txt"))
})

test_that("read_fasta_alingnment validation succeeds", {
  expect_s3_class(
    read_fasta_alignment("valid_alignment.txt"), class = "rprimer_alignment"
  )
})

test_that("remove_gaps stops when it should", {
  expect_error(remove_gaps(example_rprimer_alignment, threshold = 0.1))
  expect_error(remove_gaps(example_rprimer_alignment, threshold = 1))
  expect_error(remove_gaps(unclass(example_rprimer_alignment)))
})

test_that("remove_gaps validation succeeds", {
  expect_s3_class(
    remove_gaps(example_rprimer_alignment, threshold = 0.5),
    class = "rprimer_alignment"
  )
})

test_that("select_roi stops when it should", {
  expect_error(select_roi(example_rprimer_alignment, from = -1))
  expect_error(select_roi(example_rprimer_alignment, from = 1, to = 10))
  expect_error(select_roi(example_rprimer_alignment, from = 500, to = 10))
  expect_error(select_roi(example_rprimer_oligo))
})
