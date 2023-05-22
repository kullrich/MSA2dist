data(hiv)

test_that("rcpp_pairwiseDeletionAA()", {
    h <- hiv |> cds2aa() |> as.character()
    expect_true(rcpp_pairwiseDeletionAA(aavector=h, ncores=1)$sitesUsed[1,2]
    == 91)
    expect_true(rcpp_pairwiseDeletionAA(aavector=h, ncores=1,
    symmetric=0)$sitesUsed[1,2] == 91)
})
