context("import alignment")

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
