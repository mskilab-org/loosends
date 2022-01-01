library(testthat)
library(loosends)
library(gUtils)

tbam = system.file("tests", "testthat", "data", "pipeline", "tumor.bam", package = "loosends")
nbam = system.file("tests", "testthat", "data", "pipeline", "normal.bam", package = "loosends")
ranges = parse.gr(c("12:45763007-", "20:60158837-"))

ref.dir = "~/git/loosends/inst/extdata/hg19_loosends"
outdir = "~/testing_tmp"

test_that(desc = "process loose ends", code = {
    suppressWarnings(expr = {
        res = process_loose_ends(id = "G32831.HCC1954",
                                 ranges = ranges,
                                 tbam = tbam, nbam = nbam,
                                 ref_dir = ref.dir,
                                 outdir = outdir,
                                 read_pad = 5000,
                                 assembly_pad = 1000,
                                 verbose = FALSE)
        expect_true(!is.null(res$calls))
        expect_true(!is.null(res$contigs))
        expect_true(!is.null(res$support))
        expect_true(res$calls[, .N] == 2)
        expect_true(res$contigs[, .N] > 0)
        expect_true(res$support[, .N] == length(unique(res$contigs[, qname])))
    })
})
