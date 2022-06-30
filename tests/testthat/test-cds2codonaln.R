data(hiv)

test_that("cds2codonaln() outputs DNAStringSet", {
    cds1 <- Biostrings::DNAString("ATGCAACATTGC")
    cds2 <- Biostrings::DNAString("ATGCATTGC")
    cds1.cds2.aln <- cds2codonaln(cds1, cds2)
    expect_true(as.character(cds1.cds2.aln[1]) == "ATGCAACATTGC")
    expect_true(as.character(cds1.cds2.aln[2]) == "ATG---CATTGC")
})
