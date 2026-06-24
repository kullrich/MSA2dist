data(iupac)

test_that("dnastring2weightedPi()", {
    w1 <- dnastring2weightedPi(iupac, model="IUPAC")
    w2 <- dnastring2weightedPi(iupac, model="IUPAC", estimator="prob")
    w3 <- dnastring2weightedPi(iupac, model="sequence")
    w4 <- dnastring2weightedPi(iupac, model="sequence", estimator="prob")
    poplist <- list(FRA = grep("Mmd.FRA", names(iupac)),
        GER = grep("Mmd.GER", names(iupac)),
        IRA = grep("Mmd.IRA", names(iupac)),
        AFG = grep("Mmm.AFG", names(iupac)))
    iupac <- addpop2string(iupac, poplist)
    w5 <- dnastring2weightedPi(iupac, pop=popinteger(iupac), model="IUPAC")
    w6 <- dnastring2weightedPi(iupac, pop=popnames(iupac), model="IUPAC")
    w7 <- dnastring2weightedPi(iupac, pop=popinteger(iupac), model="IUPAC",
        estimator="prob")
    w8 <- dnastring2weightedPi(iupac, pop=popinteger(iupac), model="sequence")
    w9 <- dnastring2weightedPi(iupac, pop=popinteger(iupac), model="sequence",
        estimator="prob")
    mask1 <- IRanges::IRanges(start=c(1,61,121), end=c(30,90,150))
    w10 <- dnastring2weightedPi(iupac, model="IUPAC", mask=mask1)
    region1 <- IRanges::IRanges(start=c(1,139), end=c(75,225))
    w11 <- dnastring2weightedPi(iupac, model="IUPAC", region=region1)
    w12 <- dnastring2weightedPi(iupac, model="IUPAC", mask=mask1, region=region1)
    expect_true(signif(w1$weightedPi[[1]]) == 0.00704237)
    expect_true(signif(w2$weightedPi[[1]]) == 0.006925)
    expect_true(signif(w3$weightedPi[[1]]) == 0.00454253)
    expect_true(signif(w4$weightedPi[[1]]) == 0.00588341)
    expect_true(signif(w5$weightedPi[[1]]) == 0.00245833)
    expect_true(signif(w6$weightedPi[[1]]) == 0.00245833)
    expect_true(signif(w7$weightedPi[[1]]) == 0.00230469)
    expect_true(signif(w8$weightedPi[[1]]) == 0.00233333)
    expect_true(signif(w9$weightedPi[[1]]) == 0.00194444)
    expect_true(signif(w10$weightedPi[[1]]) == 0.00718011)
    expect_true(signif(w11$weightedPi[[1]]) == 0.0150276)
    expect_true(signif(w12$weightedPi[[1]]) == 0.0183427)
})
