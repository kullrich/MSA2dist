data(iupac)

test_that("collapseXStringSet()", {
    iupac_collapsed <- iupac |> collapseXStringSet()
    expect_true(length(iupac_collapsed$unique) == 21)
})
