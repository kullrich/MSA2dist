data(iupac)

test_that("addpos2string() outputs DNAStringSet", {
    iupac <- iupac |> addpos2string(chrom="chr1", start=1, end=1000)
    expect_true(as.character(iupac@metadata$GRanges@seqnames) == "chr1")
})
