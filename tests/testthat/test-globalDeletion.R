test_that("globalDeletion()", {
    cds1 <- Biostrings::DNAString("ATGCAACATTGC")
    cds2 <- Biostrings::DNAString("ATG---CATTGC")
    cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
        Biostrings::DNAStringSet(cds2))
    expect_true(as.character(globalDeletion(cds1.cds2.aln)[1]) == "ATGCATTGC")
})
