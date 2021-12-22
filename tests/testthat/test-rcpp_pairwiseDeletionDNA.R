data(woodmouse, package="ape")

test_that("rcpp_pairwiseDeletionDNA()", {
    w <- woodmouse |> dnabin2dnastring() |> as.character()
    expect_true(rcpp_pairwiseDeletionDNA(dnavector=w, ncores=1)$sitesUsed[2,2]
    == 962)
})
