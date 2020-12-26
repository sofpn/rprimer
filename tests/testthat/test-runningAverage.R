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
