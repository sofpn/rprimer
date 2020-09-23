test_that("gcContent works", {
  expect_equal(gcContent("ACTCTC"), 0.5)
  expect_equal(gcContent("ACTCTC--"), 0.5)
  expect_equal(gcContent("actctc--"), 0.5)
})

test_that("gcContent returns an error when it should", {
  expect_error(gcContent("ACGTM"))
  expect_error(gcContent(c("A", "C", "C")))
})
