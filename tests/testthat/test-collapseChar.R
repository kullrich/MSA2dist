data(iupac)

test_that("collapseChar()", {
    iupac_collapsed <- iupac |> as.character() |> collapseChar()
    expect_true(length(iupac_collapsed$unique) == 21)
})
