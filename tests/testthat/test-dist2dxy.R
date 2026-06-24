data(iupac)

test_that("dist2dxy()", {
    d <- dnastring2dist(iupac)
    dxy1 <- dist2dxy(d)
    poplist <- list(FRA = grep("Mmd.FRA", names(iupac)),
        GER = grep("Mmd.GER", names(iupac)),
        IRA = grep("Mmd.IRA", names(iupac)),
        AFG = grep("Mmm.AFG", names(iupac)))
    iupac <- iupac |> addpop2string(poplist)
    dxy2 <- dist2dxy(d, pop=popinteger(iupac))
    dxy3 <- dist2dxy(d, pop=popinteger(iupac), popx="FRA",
        popy="IRA", popout="AFG")
    expect_true(signif(dxy1$pi_within[1,"mean"]) == 0.00694919)
    expect_true(signif(dxy2$dxy["FRA", "GER"]) == 0.00395776)
    expect_true(signif(dxy3$RND) == 0.398615)
})
