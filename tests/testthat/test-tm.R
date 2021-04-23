# .nn ==========================================================================

test_that(".nn works", {
    x <- c("A", "C", "G")
    expect_equal(.nn(x), c("AC", "CG"))
})

# .stack =======================================================================

test_that(".stack works", {
    x <- t(matrix(.nn(c("A", "C", "G"))))
    expect_true(is.matrix(.stack(x)))
})

# .initiate ====================================================================

test_that(".initiate works", {
    x1 <- t(matrix(.nn(c("A", "T", "C"))))
    x2 <- t(matrix(.nn(c("T", "T", "C"))))
    x3 <- t(matrix(.nn(c("C", "T", "C"))))
    x4 <- t(matrix(.nn(c("A", "C", "A"))))
    x5 <- t(matrix(.nn(c("A", "C", "T"))))
    initiation <- lookup$nn$dH[lookup$nn$bases == "Initiation"]
    at_penalty <- lookup$nn$dH[lookup$nn$bases == "AT_penalty"]
    expect_equal(.initiate(x1), initiation + at_penalty)
    expect_equal(.initiate(x2), initiation + at_penalty)
    expect_equal(.initiate(x3), initiation)
    expect_equal(.initiate(x4), initiation + 2*at_penalty)
    expect_equal(.initiate(x5), initiation + 2*at_penalty)
})

# .tmParameters ================================================================

test_that(".tmParameters works", {
    x <- t(matrix(c("A", "G", "T", "T", "C", "G", "G", "T", "C", "G")))
    test <- .tmParameters(x)
    expect_true(is.matrix(test))
    test <- .tmParameters(rbind(x, x))
    expect_true(all(test[1, ] == test[2, ]))
})

# .tm ==========================================================================

test_that(".tm works", {
    x <- .tmParameters(t(matrix(
        c("A", "G", "T", "T", "C", "G", "G", "T", "C")
        )))
    test1 <- .tm(x, 250)
    expect_true(is.numeric(test1))
    test2 <- .tm(x, 500)
    expect_true(is.numeric(test2))
    expect_true(test2 > test1)
    x <- rbind(x, x)
    test <- .tm(x)
    expect_true(test[1] == test[2])
})

# .deltaG ======================================================================

#test_that(".deltaG works", {
#    ## Confirm that the calculation agrees with the examples in
    ## SantaLucia and Hicks 2004
#    m <- t(matrix(unlist(strsplit("CGTTGA", split = ""))))
#    x <- .tmParameters(m, concNa = 1)
#    test <- .deltaG(x)
#    expect_equal(round(test, 2), -5.36)
#    x <- .tmParameters(m, concNa = 0.115)
#    test <- .deltaG(x)
#    expect_equal(round(test, 2), -4.12)

#})
