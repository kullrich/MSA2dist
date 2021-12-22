data(iupac)

test_that("addmask2string() outputs DNAStringSet", {
    mask1 <- IRanges::IRanges(start=c(1,41), end=c(20,50))
    iupac <- iupac |> addmask2string(mask=mask1)
    expect_true(iupac@metadata$mask@width[1] == 20)
    mask2 <- IRanges::IRanges(start=c(21), end=c(30))
    iupac <- iupac |> addmask2string(mask=mask2)
    expect_true(iupac@metadata$mask@width[1] == 30)
    iupac <- iupac |> addmask2string(mask=mask2, append=FALSE)
    expect_true(iupac@metadata$mask@width[1] == 10)
})
