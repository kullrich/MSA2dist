data(woodmouse, package="ape")

test_that("aa2selfscore() outputs data.frame", {
    aa <- woodmouse |> dnabin2dnastring() |> cds2aa(shorten=TRUE,
        genetic.code=Biostrings::getGeneticCode("2"))
    expect_true(class(aa) == "AAStringSet")
    expect_true((aa |> aa2selfscore(scorematrix = "BLOSUM45"))[1, 2] == 2043)
    expect_true((aa |> aa2selfscore(scorematrix = "BLOSUM50"))[1, 2] == 2156)
    expect_true((aa |> aa2selfscore(scorematrix = "BLOSUM62"))[1, 2] == 1690)
    expect_true((aa |> aa2selfscore(scorematrix = "BLOSUM80"))[1, 2] == 1809)
    expect_true((aa |> aa2selfscore(scorematrix = "BLOSUM100"))[1, 2] == 3132)
    expect_true((aa |> aa2selfscore(scorematrix = "PAM30"))[1, 2] == 2484)
    expect_true((aa |> aa2selfscore(scorematrix = "PAM40"))[1, 2] == 2438)
    expect_true((aa |> aa2selfscore(scorematrix = "PAM70"))[1, 2] == 2178)
    expect_true((aa |> aa2selfscore(scorematrix = "PAM120"))[1, 2] == 1788)
    expect_true((aa |> aa2selfscore(scorematrix = "PAM250"))[1, 2] == 1747)
})
