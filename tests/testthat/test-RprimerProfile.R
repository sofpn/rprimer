test_that("RprimerProfile validation works", {
  expect_error(RprimerProfile(matrix(2, 1, 1)))
})
