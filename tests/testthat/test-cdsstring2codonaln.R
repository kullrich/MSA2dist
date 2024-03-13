data(hiv)

test_that("cdsstring2codonaln() outputs DNAStringSet", {
    cds <- Biostrings::DNAStringSet(c("ATGCAACATTGC", "ATGCATTGC"))
    names(cds) <- c("cds1", "cds2")
    aa <- MSA2dist::cds2aa(cds)
    cds1.cds2.aln <- cdsstring2codonaln(cds, aa)
    expect_true(as.character(cds1.cds2.aln[1]) == "ATGCAACATTGC")
    expect_true(as.character(cds1.cds2.aln[2]) == "ATG---CATTGC")
})
