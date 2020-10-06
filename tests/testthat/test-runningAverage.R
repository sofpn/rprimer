test_that(".runningAverage works", {
    toTest <- runif(100)
    expect_identical(nrow(.runningAverage(toTest)), length(toTest))
})

test_that(".runningAverage returns an error when it should", {
    toTest <- runif(100)
    expect_error(.runningAverage(toTest, size = 300))
    expect_error(.runningAverage(toTest, size = 0))
})

test_that(".gcRunningAverage works", {
    validSequence <- c("A", "C", "T", "G", "T", "T", "G", "C", "A")
    res <- .gcRunningAverage(validSequence)
    expect_equal(nrow(res), length(validSequence))
})

test_that(".gcRunningAverage returns an error when it should", {
    invalidSequence <- c("A", "C", "T", "G", "M", "T", "G", "C", "A")
    validSequence <- c("A", "C", "T", "G", "T", "T", "G", "C", "A")
    expect_error(.gcRunningAverage(invalidSequence))
    expect_error(.gcRunningAverage(validSequence, size = 300))
    expect_error(.gcRunningAverage(validSequence, size = 0))
})
