data(hiv)

test_that("rcpp_KaKs()", {
    h <- rcpp_KaKs(cdsstr=as.character(hiv[1:3]))
    expect_true(h[[3]][1] == "U68496_U68497")
})
