data(iupac)

test_that("dnastring2dxy()", {
    w1 <- dnastring2dxy(iupac, model="IUPAC")
    w2 <- dnastring2dxy(iupac, model="IUPAC", estimator="count")
    w3 <- dnastring2dxy(iupac, model="sequence")
    w4 <- dnastring2dxy(iupac, model="sequence", estimator="count")
    poplist <- list(FRA = grep("Mmd.FRA", names(iupac)),
        GER = grep("Mmd.GER", names(iupac)),
        IRA = grep("Mmd.IRA", names(iupac)),
        AFG = grep("Mmm.AFG", names(iupac)))
    iupac <- addpop2string(iupac, poplist)
    w5 <- dnastring2dxy(iupac, pop=popinteger(iupac), model="IUPAC")
    w6 <- dnastring2dxy(iupac, pop=popnames(iupac), model="IUPAC")
    w7 <- dnastring2dxy(iupac, pop=popinteger(iupac), model="IUPAC",
        estimator="count")
    w8 <- dnastring2dxy(iupac, pop=popinteger(iupac), model="sequence")
    w9 <- dnastring2dxy(iupac, pop=popinteger(iupac), model="sequence",
        estimator="count")
    mask1 <- IRanges::IRanges(start=c(1,61,121), end=c(30,90,150))
    w10 <- dnastring2dxy(iupac, model="IUPAC", mask=mask1)
    region1 <- IRanges::IRanges(start=c(1,139), end=c(75,225))
    w11 <- dnastring2dxy(iupac, model="IUPAC", region=region1)
    w12 <- dnastring2dxy(iupac, model="IUPAC", mask=mask1, region=region1)
    expect_true(signif(w1$dxy[[1]]) == 0.006925)
    expect_true(signif(w2$dxy[[1]]) == 0.00704237)
    expect_true(signif(w3$dxy[[1]]) == 0.00588341)
    expect_true(signif(w4$dxy[[1]]) == 0.00612497)
    expect_true(signif(w5$dxy[[1]]) == 0.00230469)
    expect_true(signif(w6$dxy[[1]]) == 0.00230469)
    expect_true(signif(w7$dxy[[1]]) == 0.00245833)
    expect_true(signif(w8$dxy[[1]]) == 0.00194444)
    expect_true(signif(w9$dxy[[1]]) == 0.00233333)
    expect_true(signif(w10$dxy[[1]]) == 0.00706044)
    expect_true(signif(w11$dxy[[1]]) == 0.0147771)
    expect_true(signif(w12$dxy[[1]]) == 0.018037)
})
