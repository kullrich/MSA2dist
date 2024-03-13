data(hiv)

test_that("compareCodons()", {
    idx <- list(c(1,2), c(1,3))
    h <- hiv[1:3] |> indices2kaks(idx, model="Li")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="Li",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="NG86")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="NG86",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="NG")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="NG",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="LWL")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="LWL",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="LPB")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="LPB",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="MLWL")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="MLWL",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="MLPB")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="MLPB",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="GY")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="GY",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="YN")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="YN",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="MYN")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="MYN",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="MS")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="MS",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="MA")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="MA",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="GNG")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="GNG",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="GLWL")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="GLWL",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="GLPB")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="GLPB",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="GMLWL")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="GMLWL",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="GMLPB")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="GMLPB",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="GYN")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="GYN",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="GMYN")
    expect_true(h$seq1[1] == "U68496")
    h <- hiv[1:3] |> indices2kaks(idx, model="GMYN",
        isMSA=FALSE)
    expect_true(h$seq1[1] == "U68496")
})
