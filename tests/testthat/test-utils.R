## library(gUtils)
## library(loosends)
## library(testthat)

## human.fasta.fn = system.file("extdata", "hg19_loosends", "human_g1k_v37_decoy.fasta", package = "loosends")

## test_that(desc = "check that minimap pre-indexing works", code = {
##     suppressWarnings(
##         expr = {
##             mmi = minimap_index(fasta = human.fasta.fn,
##                                 outdir = "~/testing_tmp/index",
##                                 verbose = FALSE)
##             ## make sure that a non-empty file is created
##             expect_true(file.exists(mmi))
##             expect_true(file.info(mmi)$size > 0)
##         }
##     )
## })

