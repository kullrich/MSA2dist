data(iupac)

test_that("popinteger() outputs list", {
    poplist <- list(FRA = grep("Mmd.FRA", names(iupac)))
    iupac <- iupac |> addpop2string(poplist)
    expect_true(popinteger(iupac)$FRA[1] == 1)
})
