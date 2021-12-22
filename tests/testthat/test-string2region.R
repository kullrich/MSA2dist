data(iupac)

test_that("string2region() outputs DNAStringSet", {
    mask1 <- IRanges::IRanges(start=c(11,41,71), end=c(20,50,80))
    hiv.region <- hiv |> cds2aa() |> string2region(mask=mask1)
    expect_true(regionused(hiv.region)@start[1] == 1)
})
