data(hiv)

test_that("rcpp_distSTRING()", {
    expect_true(rcpp_distSTRING(dnavector=as.character(hiv),
    scoreMatrix=iupacMatrix())$sitesUsed[1,2] == 273)
    expect_true(rcpp_distSTRING(dnavector=as.character(hiv),
    scoreMatrix=iupacMatrix(), ncores=2)$sitesUsed[1,2] == 273)
    myscore <- iupacMatrix()
    myscore[1,4] <- 0.5
    h <- rcpp_distSTRING(dnavector=as.character(hiv), scoreMatrix=myscore)
    expect_true(h$distSTRING[1,2] == h$distSTRING[2,1])
    h <- rcpp_distSTRING(dnavector=as.character(hiv), scoreMatrix=myscore,
        symmetric=0)
    expect_true(h$distSTRING[1,2] != h$distSTRING[2,1])
})
