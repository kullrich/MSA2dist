data(iupac)

test_that("addregion2string() outputs DNAStringSet", {
    region1 <- IRanges::IRanges(start=c(1,41), end=c(20,50))
    iupac <- iupac |> addregion2string(region=region1)
    expect_true(iupac@metadata$region@width[1] == 20)
    region2 <- IRanges::IRanges(start=c(21), end=c(30))
    iupac <- iupac |> addregion2string(region=region2)
    expect_true(iupac@metadata$region@width[1] == 30)
    iupac <- iupac |> addregion2string(region=region2, append=FALSE)
    expect_true(iupac@metadata$region@width[1] == 10)
})
