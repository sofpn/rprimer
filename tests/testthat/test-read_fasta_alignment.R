context("import alignment")

test_that("read_fasta_alignment stops when it should", {
  expect_error(read_fasta_alignment(invalid_filename))
  expect_error(read_fasta_alignment("file_not_found.txt"))
  expect_error(read_fasta_alignment("non_fasta_format.txt"))
  expect_error(read_fasta_alignment("not_unique_sequence_names.txt"))
  expect_error(read_fasta_alignment("different_sequence_lengths.txt"))
  expect_error(read_fasta_alignment("invalid_sequence_character.txt"))
})

test_that("read_fasta_alingnment validation succeeds", {
  expect_s3_class(
    read_fasta_alignment("valid_alignment.txt"),
    class = "rprimer_alignment"
  )
})
