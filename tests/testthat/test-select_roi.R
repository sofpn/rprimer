context("import alignment")

test_that("select_roi stops when it should", {
  expect_error(select_roi(example_rprimer_alignment, from = -1))
  expect_error(select_roi(example_rprimer_alignment, from = 500, to = 10))
  expect_error(select_roi(example_rprimer_oligo))
})
