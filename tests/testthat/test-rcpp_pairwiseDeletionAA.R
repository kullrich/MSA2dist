data(hiv)

test_that("rcpp_pairwiseDeletionAA()", {
    h <- hiv |> cds2aa() |> as.character()
    expect_true(rcpp_pairwiseDeletionAA(aavector=h, ncores=1)$sitesUsed[1,2]
    == 91)
})
