data(hiv)

test_that("rcpp_weightedPi()", {
    expect_true(signif(rcpp_weightedPi(dnavector=as.character(hiv),
    model="IUPAC", estimator="count")$weightedPi) == 0.096523)
    expect_true(signif(rcpp_weightedPi(dnavector=as.character(hiv),
    model="IUPAC", estimator="prob")$weightedPi) == 0.0928105)
    expect_true(signif(rcpp_weightedPi(dnavector=as.character(hiv),
    model="sequence", estimator="count")$weightedPi) == 0.100545)
    expect_true(signif(rcpp_weightedPi(dnavector=as.character(hiv),
    model="sequence", estimator="prob")$weightedPi) == 0.0928105)
})
