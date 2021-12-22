data(woodmouse, package="ape")

test_that("aabin2aastring() outputs AAStringSet", {
    expect_true(class(ape::trans(woodmouse[,1:6], 2) |> aabin2aastring())
    == "AAStringSet")
    expect_true(class(Biostrings::AAStringSet("M") |> aastring2aabin() |>
    aabin2aastring()) == "AAStringSet")
})
