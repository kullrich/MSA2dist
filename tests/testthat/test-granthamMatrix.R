test_that("granthamMatrix()", {
    expect_true(granthamMatrix()[1,10] == 155)
    expect_true(granthamMatrix()[6,7] == 64)
})
