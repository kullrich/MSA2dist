data(iupac)

test_that("rcpp_weightedPi_pop()", {
    weightedPi_1 <- rcpp_weightedPi_pop(dnavector=as.character(iupac),
    pop_idx=c(rep(1,8),rep(2,8),rep(3,8),rep(4,6)),
    model="IUPAC",
    estimator="count")
    expect_true(signif(weightedPi_1$weightedPi[[1]]) == 0.00245833)
    expect_true(signif(weightedPi_1$weightedPi[[2]]) == 0.00409167)
    expect_true(signif(weightedPi_1$weightedPi[[3]]) == 0.00674167)
    expect_true(signif(weightedPi_1$weightedPi[[4]]) == 0.000941276)
    weightedPi_2 <- rcpp_weightedPi_pop(dnavector=as.character(iupac),
    pop_idx=c(rep(1,8),rep(2,8),rep(3,8),rep(4,6)),
    model="IUPAC",
    estimator="prob")
    expect_true(signif(weightedPi_2$weightedPi[[1]]) == 0.00230469)
    expect_true(signif(weightedPi_2$weightedPi[[2]]) == 0.00383594)
    expect_true(signif(weightedPi_2$weightedPi[[3]]) == 0.00632031)
    expect_true(signif(weightedPi_2$weightedPi[[4]]) == 0.000862837)
    weightedPi_3 <- rcpp_weightedPi_pop(dnavector=as.character(iupac),
    pop_idx=c(rep(1,8),rep(2,8),rep(3,8),rep(4,6)),
    model="sequence",
    estimator="count")
    expect_true(signif(weightedPi_3$weightedPi[[1]]) == 0.00233333)
    expect_true(signif(weightedPi_3$weightedPi[[2]]) == 0.00411905)
    expect_true(signif(weightedPi_3$weightedPi[[3]]) == 0.0044)
    expect_true(signif(weightedPi_3$weightedPi[[4]]) == 0)
    weightedPi_4 <- rcpp_weightedPi_pop(dnavector=as.character(iupac),
    pop_idx=c(rep(1,8),rep(2,8),rep(3,8),rep(4,6)),
    model="sequence",
    estimator="prob")
    expect_true(signif(weightedPi_4$weightedPi[[1]]) == 0.00194444)
    expect_true(signif(weightedPi_4$weightedPi[[2]]) == 0.00357526)
    expect_true(signif(weightedPi_4$weightedPi[[3]]) == 0.00304222)
    expect_true(signif(weightedPi_4$weightedPi[[4]]) == 0)
})
