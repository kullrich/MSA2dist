test_that("uptriidx()", {
    expect_true(uptriidx(3)[3] == 6)
    expect_true(uptriidx(2, diag=TRUE)[3] == 4)
})
