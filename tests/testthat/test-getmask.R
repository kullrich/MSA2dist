data(iupac)

test_that("getmask() outputs IRanges object", {
    mask1 <- IRanges::IRanges(start=c(1,41), end=c(20,50))
    iupac <- iupac |> addmask2string(mask=mask1)
    expect_true(getmask(iupac)@width[1] == 20)
})
