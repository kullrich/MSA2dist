data(iupac)

test_that("getpos() outputs GRanges object", {
    iupac <- iupac |> addpos2string(chrom="chr1", start=1, end=1000)
    expect_true(as.character(getpos(iupac)@seqnames) == "chr1")
})
