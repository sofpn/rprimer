test_that("runningAverage works", {
  toTest <- runif(100)
  expect_identical(nrow(runningAverage(toTest)), length(toTest))
})

test_that("runningAverage returns an error when it should", {
  toTest <- runif(100)
  expect_error(runningAverage(toTest, size = 300))
  expect_error(runningAverage(toTest, size = 0))
})
