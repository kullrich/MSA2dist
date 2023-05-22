data(hiv)

test_that("aastring2dist()", {
    h <- hiv |> cds2aa() |> aastring2dist(score=granthamMatrix())
    expect_true(h$sitesUsed[1,1] == 91)
    h <- hiv |> cds2aa() |> aastring2dist(score=granthamMatrix(), threads=2)
    expect_true(h$sitesUsed[1,1] == 91)
    mask1 <- IRanges::IRanges(start=c(11,41,71), end=c(20,50,80))
    h <- hiv |> cds2aa() |> aastring2dist(score=granthamMatrix(), mask=mask1)
    expect_true(h$sitesUsed[1,1] == 61)
    region1 <- IRanges::IRanges(start=c(1,75), end=c(45,85))
    h <- hiv |> cds2aa() |> aastring2dist(score=granthamMatrix(),
        region=region1)
    expect_true(h$sitesUsed[1,1] == 56)
    h <- hiv |> cds2aa() |> aastring2dist(score=granthamMatrix(),
        mask=mask1, region=region1)
    expect_true(h$sitesUsed[1,1] == 35)
    myscore <- granthamMatrix()
    myscore[5,6]<-0
    h <- hiv |> cds2aa() |> aastring2dist(score=myscore)
    expect_true(h$distSTRING[1,2] == h$distSTRING[2,1])
    h <- hiv |> cds2aa() |> aastring2dist(score=myscore, symmetric=FALSE)
    expect_true(h$distSTRING[1,2] != h$distSTRING[2,1])
})
