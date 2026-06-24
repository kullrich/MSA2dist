data(iupac)

test_that("rcpp_dxy_fst_pop()", {
    dxy_1 <- rcpp_dxy_fst_pop(dnavector=as.character(iupac),
    pop_idx=c(rep(1,8),rep(2,8),rep(3,8),rep(4,6)),
    model="IUPAC",
    estimator="count")
    expect_true(signif(dxy_1$dxy[1,1]) == 0.00245833)
    expect_true(signif(dxy_1$dxy[1,2]) == 0.00401562)
    expect_true(signif(dxy_1$dxy[1,3]) == 0.00517188)
    expect_true(signif(dxy_1$dxy[1,4]) == 0.0122537)
    dxy_2 <- rcpp_dxy_fst_pop(dnavector=as.character(iupac),
    pop_idx=c(rep(1,8),rep(2,8),rep(3,8),rep(4,6)),
    model="IUPAC",
    estimator="prob")
    expect_true(signif(dxy_2$dxy[1,1]) == 0.00230469)
    expect_true(signif(dxy_2$dxy[1,2]) == 0.00401562)
    expect_true(signif(dxy_2$dxy[1,3]) == 0.00517188)
    expect_true(signif(dxy_2$dxy[1,4]) == 0.0122537)
    dxy_3 <- rcpp_dxy_fst_pop(dnavector=as.character(iupac),
    pop_idx=c(rep(1,8),rep(2,8),rep(3,8),rep(4,6)),
    model="sequence",
    estimator="count")
    expect_true(signif(dxy_3$dxy[1,1]) == 0.00233333)
    expect_true(signif(dxy_3$dxy[1,2]) == 0.00384524)
    expect_true(signif(dxy_3$dxy[1,3]) == 0.00326667)
    expect_true(signif(dxy_3$dxy[1,4]) == 0.012191)
    dxy_4 <- rcpp_dxy_fst_pop(dnavector=as.character(iupac),
    pop_idx=c(rep(1,8),rep(2,8),rep(3,8),rep(4,6)),
    model="sequence",
    estimator="prob")
    expect_true(signif(dxy_4$dxy[1,1]) == 0.00194444)
    expect_true(signif(dxy_4$dxy[1,2]) == 0.00384524)
    expect_true(signif(dxy_4$dxy[1,3]) == 0.00326667)
    expect_true(signif(dxy_4$dxy[1,4]) == 0.012191)
})
