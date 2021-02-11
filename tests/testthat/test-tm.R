test_that(".nn works", {
    x <- c("A", "C", "G")
    expect_equal(.nn(x), c("AC", "CG"))
})

test_that(".stack works", {
    x <- t(matrix(.nn(c("A", "C", "G"))))
    expect_true(is.matrix(.stack(x)))
})
