data(hiv)

test_that("rcpp_distSTRING()", {
    expect_true(rcpp_distSTRING(dnavector=as.character(hiv),
    scoreMatrix=iupacMatrix())$sitesUsed[1,2] == 273)
    expect_true(rcpp_distSTRING(dnavector=as.character(hiv),
    scoreMatrix=iupacMatrix(), ncores=2)$sitesUsed[1,2] == 273)
})
