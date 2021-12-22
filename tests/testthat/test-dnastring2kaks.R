data(hiv)

test_that("compareCodons()", {
    h <- hiv |> dnastring2kaks(model="Li")
    expect_true(h$Comp1[1] == "U68496")
    h <- hiv |> dnastring2kaks(model="NG86")
    expect_true(h$Comp1[1] == "1")
    h <- hiv |> dnastring2kaks(model="NG86", threads=2)
    expect_true(h$seq1[1] == "U68496")
})
