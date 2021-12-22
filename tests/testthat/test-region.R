data(iupac)

test_that("region() outputs IRanges object", {
    region1 <- IRanges::IRanges(start=c(1,41), end=c(20,50))
    iupac <- iupac |> addregion2string(region=region1)
    expect_true(region(iupac)@width[1] == 20)
})
