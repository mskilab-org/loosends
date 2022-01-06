library(testthat)
library(loosends)

tbam = system.file("extdata", "tests", "pipeline", "tumor.bam", package = "loosends")
nbam = system.file("extdata", "tests", "pipeline", "normal.bam", package = "loosends")
ranges = parse.gr(c("12:45763007-", "20:60158837-"))

ref.dir = "~/git/loosends/inst/extdata/hg19_loosends"
outdir = "~/testing_tmp"

res = suppressWarnings(
    process_loose_ends(id = "G32831.HCC1954",
                       ranges = ranges,
                       tbam = tbam, nbam = nbam,
                       ref_dir = ref.dir,
                       outdir = outdir,
                       read_pad = 5000,
                       assembly_pad = 1000,
                       verbose = FALSE)
)

test_that(desc = "process loose ends", code = {
    suppressWarnings(expr = {
        expect_true(!is.null(res$calls))
        expect_true(!is.null(res$contigs))
        expect_true(!is.null(res$support))
        expect_true(res$calls[, .N] == 2)
        expect_true(res$contigs[, .N] > 0)
        expect_true(res$support[, .N] == length(unique(res$contigs[, qname])))
    })
})
