test_that("countDegenerates works", {
  expect_equal(countDegenerates("CRTN"), 2)
  expect_equal(countDegenerates("crtn"), 2)
  expect_equal(countDegenerates("CRTN-"), 2)
})

test_that("countDegenerates returns an error when it should", {
  expect_error(countDegenerates(c("A","C")))
  expect_error(countDegenerates("X"))

})
