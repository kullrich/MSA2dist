data(hiv)

test_that("codonmat2pnps()", {
    expect_true(((hiv |> dnastring2codonmat())[,c(1, 2)] |>
    codonmat2pnps())["Codons"] == "91")
})
