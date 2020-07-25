context("get oligos")

test_that("get_nmers returns an error when it should", {
  expect_error(get_nmers(c(1, 2, 4)))
  expect_error(get_nmers(c("g", "t", "a", "g", "c", "g"), n = 20))
  expect_error(get_nmers(c("g", "t", "a", "g", "c", "g"), n = 0.5))
})

test_that("get_nmers works", {
  vec <- c("g", "t", "r", "a", "g", "c")
  expect_equal(get_nmers(vec), c("g", "t", "r", "a", "g", "c"))
  expect_equal(get_nmers(vec, n = 2), c("gt", "tr", "ra", "ag", "gc"))
  expect_equal(get_nmers(vec, n = 6), "gtragc")
})

test_that("gc_content returns an error when it should", {
  expect_error(gc_content("crgga"))
  expect_error(gc_content("c", "g", "g", "a"))
  expect_error(gc_content(2))
})

test_that("gc_content works", {
  expect_equal(gc_content("cgtt"), 0.5)
  expect_equal(gc_content("cgtt--"), 0.5)
  expect_equal(gc_content("CGTT"), 0.5)
})

test_that("complement returns an error when it should", {
  expect_error(complement(2))
  expect_error(complement(c("c", "g", "t")))
  expect_error(complement("cgtrzx"))
})

test_that("complement works", {
  expect_equal(complement("c"), "g")
  expect_equal(complement("ctt"), "gaa")
  expect_equal(complement("r"), "y")
  expect_equal(complement("R"), "y")
})

test_that("reverse_complement returns an error when it should", {
  expect_error(reverse_complement(2))
  expect_error(reverse_complement(c("g", "t")))
  expect_error(reverse_complement("gxttn"))
})

test_that("reverse_complement works", {
  expect_equal(reverse_complement("acg"), "cgt")
  expect_equal(reverse_complement("ACG"), "cgt")
  expect_equal(reverse_complement("yrn"), "nyr")
  expect_equal(reverse_complement("-a"), "t-")
})

test_that("running_sum returns an error when it should", {
  expect_error(running_sum("c"))
  expect_error(running_sum(example_rprimer_sequence_profile))
  expect_error(running_sum(c(2, 3, 5, 4), n = 0))
  expect_error(running_sum(c(2, 3, 5, 4), n = 10))
})

test_that("running_sum works", {
  set.seed(3)
  x <- runif(20)
  expect_equal(running_sum(x, n = 3)[[1]], x[[1]] + x[[2]] + x[[3]])
  expect_true(is.double(running_sum(x)))
  expect_equal(length(running_sum(x)), length(x) - 1)
  expect_equal(length(running_sum(x, n = 3)), length(x) - 2)
})

test_that("repeat_pattern returns an error when it should", {
  expect_error(repeat_pattern(0))
  expect_error(repeat_pattern("a"))
})

test_that("repeat_pattern works", {
  expect_equal(repeat_pattern(3), "\\1\\1")
})


test_that("exclude_oligos returns an error when it should", {
 expect_error(exclude_oligos(12))
 expect_error(exclude_oligos("cgg", pattern = 12))
})

test_that("exclude_oligos works", {
 expect_true(is.na(exclude_oligos("at", "t$")))
  expect_true(is.na(exclude_oligos("attt", paste0("(t)", repeat_pattern(3)))))
})
