context("sequence profile")

test_that("sequence_profile stops when it should", {
 expect_error(sequence_profile(unclass(example_rprimer_alignment)))
})

test_that("sequence_profile works as it should", {
  res <- sequence_profile(example_rprimer_alignment)
  expect_equal(
    dim(res)[[2]],
    nchar(example_rprimer_alignment[[1]])
  )
  expect_s3_class(res, "rprimer_sequence_profile")
  expect_true(is.numeric(res))
  expect_false(any(res > 1))
  expect_false(any(res < 0))
})
