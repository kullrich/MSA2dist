data(woodmouse, package="ape")

test_that("dnabin2dnastring outputs DNAStringSet", {
    expect_true(class(woodmouse |> dnabin2dnastring())
    == "DNAStringSet")
})
