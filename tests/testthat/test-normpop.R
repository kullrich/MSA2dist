data(iupac)

test_that("normpop()", {
    poplist <- list(FRA = grep("Mmd.FRA", names(iupac)),
        GER = grep("Mmd.GER", names(iupac)),
        IRA = grep("Mmd.IRA", names(iupac)),
        AFG = grep("Mmm.AFG", names(iupac)))
    npop1 <- normpop(pop=poplist, seqnames=names(iupac))
    poplist <- list(FRA = names(iupac)[grep("Mmd.FRA", names(iupac))],
        GER = grep("Mmd.GER", names(iupac)),
        IRA = names(iupac)[grep("Mmd.IRA", names(iupac))],
        AFG = grep("Mmm.AFG", names(iupac)))
    npop2 <- normpop(pop=poplist, seqnames=names(iupac))
    poplist <- list(FRA = names(iupac)[grep("Mmd.FRA", names(iupac))],
        GER = grep("Mmd.GER", names(iupac)),
        IRA = names(iupac)[grep("Mmd.IRA", names(iupac))])
    npop3 <- normpop(pop=poplist, seqnames=names(iupac))
    popchar <- c(rep("FRA", 8),
        rep("GER", 8),
        rep("IRA", 8),
        rep("AFG", 6))
    npop4 <- normpop(pop=popchar, seqnames=names(iupac))
    popcharn <- c(rep("FRA", 8),
        rep("GER", 8),
        rep("IRA", 8),
        rep("AFG", 6))
    names(popcharn) <- names(iupac)
    npop5 <- normpop(pop=popcharn, seqnames=names(iupac))
    popint <- c(rep(1L, 8),
        rep(2L, 8),
        rep(3L, 8),
        rep(4L, 6))
    npop6 <- normpop(pop=popint, seqnames=names(iupac))
    popintn <- c(rep(1L, 8),
        rep(2L, 8),
        rep(3L, 8),
        rep(4L, 6))
    names(popintn) <- names(iupac)
    npop7 <- normpop(pop=popintn, seqnames=names(iupac))
    expect_true(npop1$pop_idx[1] == 1)
    expect_true(npop2$pop_idx[1] == 1)
    expect_true(npop3$pop_idx[1] == 1)
    expect_true(npop4$pop_idx[1] == 1)
    expect_true(npop5$pop_idx[1] == 1)
    expect_true(npop6$pop_idx[1] == 1)
    expect_true(npop7$pop_idx[1] == 1)
})
