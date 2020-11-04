data("exampleRprimerProfile")

test_that(".runningAverage works", {
    toTest <- runif(100)
    expect_identical(nrow(.runningAverage(toTest)), length(toTest))
    toTest <- runif(10)
    expect_identical(nrow(.runningAverage(toTest)), length(toTest))
    toTest <- runif(1000)
    expect_equal(.runningAverage(toTest, size = 100)$position[1], 50)
})

test_that(".runningAverage returns an error when it should", {
    toTest <- runif(100)
    expect_error(.runningAverage(toTest, size = 300))
    expect_error(.runningAverage(toTest, size = 0))
    expect_error(.runningAverage(toTest, size = FALSE))
})

test_that(".gcRunningAverage works", {
    validSequence <- c("A", "C", "T", "G", "T", "T", "G", "C", "A")
    res <- .gcRunningAverage(validSequence)
    expect_equal(nrow(res), length(validSequence))
    toTest <- exampleRprimerProfile$majority
    expect_equal(.gcRunningAverage(toTest, size = 200)$position[1], 100)
})

test_that(".gcRunningAverage returns an error when it should", {
    validSequence <- c("A", "C", "T", "G", "T", "T", "G", "C", "A")
    expect_error(.gcRunningAverage(validSequence, size = 300))
    expect_error(.gcRunningAverage(validSequence, size = 0))
    expect_error(.gcRunningAverage(validSequence, size = FALSE))
})
