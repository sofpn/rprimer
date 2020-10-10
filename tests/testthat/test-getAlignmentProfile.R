test_that("getAlignmentProfile works", {
  data("exampleRprimerAlignment")
  toTest <- getAlignmentProfile(exampleRprimerAlignment)
  expect_true(is(toTest, "RprimerProfile"))
})
