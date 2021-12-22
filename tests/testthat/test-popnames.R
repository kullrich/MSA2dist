data(iupac)

test_that("popnames() outputs list", {
    poplist <- list(FRA = grep("Mmd.FRA", names(iupac)))
    iupac <- iupac |> addpop2string(poplist)
    expect_true(popnames(iupac)$FRA[1] == "Mmd.FRA.14")
})
