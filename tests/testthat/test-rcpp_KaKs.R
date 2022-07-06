data(hiv)

test_that("rcpp_KaKs()", {
    h <- rcpp_KaKs(cdsstr=as.character(hiv[1:3]))
    expect_true(h[3,2] == "U68496")
})
