data(hiv)

test_that("compareCodons()", {
    h <- hiv[1:3] |> dnastring2kaks(model="Li")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="Li",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="NG86")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="NG86",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="NG")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="NG",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="LWL")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="LWL",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="LPB")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="LPB",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="MLWL")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="MLWL",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="MLPB")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="MLPB",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="GY")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="GY",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="YN")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="YN",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="MYN")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="MYN",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="MS")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="MS",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="MA")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="MA",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="GNG")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="GNG",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="GLWL")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="GLWL",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="GLPB")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="GLPB",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="GMLWL")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="GMLWL",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="GMLPB")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="GMLPB",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="GYN")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="GYN",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="GMYN")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> dnastring2kaks(model="GMYN",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
})
