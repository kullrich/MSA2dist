test_that("aln2dnastring() outputs DNAStringSet", {
    cds1 <- Biostrings::DNAString("ATGCAACATTGC")
    cds2 <- Biostrings::DNAString("ATG---CATTGC")
    cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
        Biostrings::DNAStringSet(cds2))
    expect_true(class(cds1.cds2.aln |> dnastring2aln() |> aln2dnastring())
    == "DNAStringSet")
})
