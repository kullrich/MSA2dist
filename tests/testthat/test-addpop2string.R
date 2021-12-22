data(iupac)

test_that("addpop2string() outputs DNAStringSet", {
    poplist <- list(FRA = grep("Mmd.FRA", names(iupac)))
    iupac <- iupac |> addpop2string(poplist)
    expect_true(iupac@metadata$pop.integer$FRA[1] == 1)
})
