data(hiv)

test_that("codonmat2pnps()", {
    expect_true((hiv |> dnastring2codonmat() |> codonmat2xy())$NonSynSum[1]
    == 40)
    expect_true((hiv |> dnastring2codonmat() |>
    codonmat2xy(threads=2))$NonSynSum[1] == 40)
})
