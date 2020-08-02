context("import alignment")

test_that("remove_gaps stops when it should", {
  expect_error(remove_gaps(example_rprimer_alignment, threshold = 0.1))
  expect_error(remove_gaps(example_rprimer_alignment, threshold = 1))
  expect_error(remove_gaps(unclass(example_rprimer_alignment)))
})

test_that("remove_gaps succeeds", {
  expect_s3_class(
    remove_gaps(example_rprimer_alignment, threshold = 0.5),
    class = "rprimer_alignment"
  )
})
