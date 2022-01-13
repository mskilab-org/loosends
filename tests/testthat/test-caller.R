library(gUtils)
library(loosends)
library(testthat)
library(RSeqLib)

## reference directory
ref.dir = "~/git/loosends/inst/extdata/hg19_loosends"
loosereads.fn = system.file("extdata", "tests", "loosereads_1", "loosereads.res.rds", package = "loosends")
reads.dt.fn = system.file("extdata", "tests", "loosereads_1", "reads.dt.rds", package = "loosends")
le.dt.fn = system.file("extdata", "tests", "loosereads_1", "le.dt.rds", package = "loosends")
big.reads.dt.fn = system.file("extdata", "tests", "loosereads_1", "big.reads.dt.rds", package = "loosends")
big.le.dt.fn = system.file("extdata", "tests", "loosereads_1", "big.le.dt.rds", package = "loosends")
concat.fasta.fn = system.file("extdata", "hg19_loosends", "concatenated_references_deduped.fasta", package = "loosends")
human.fasta.fn = system.file("extdata", "hg19_loosends", "human_g1k_v37_decoy.fasta", package = "loosends")

## useful params
this.pair = "G32831.HCC1954"
this.le = parse.gr("20:60158837-")

## read input files
loosereads.dt = readRDS(loosereads.fn)
le.dt = readRDS(le.dt.fn)
reads.dt = readRDS(reads.dt.fn)

## read input files for pipeline
big.le.dt = readRDS(big.le.dt.fn)
big.reads.dt = readRDS(big.reads.dt.fn)

## load ref bwa
## ref.obj = grab_ref_obj(ref.dir = ref.dir)
concat = BWA(fasta = concat.fasta.fn)
human = BWA(fasta = human.fasta.fn)

## test_that(desc = "check that loose ends are put into correct format", code = {
##     suppressWarnings(
##         expr = {
##             le.dt = prep_loose_ends(this.le, id = this.pair)
##             expect_true(!is.null(le.dt$leix))
##             expect_true(all(le.dt$sample == this.pair))
##         }
##     )
## })

## test_that(desc = "check that loose reads are put into the correct format", code = {
##     suppressWarnings(
##         expr = {
##             reads.dt = prep_loose_reads(li = le.dt, loose.reads.dt = loosereads.dt)
##             check that concord and anchor are added
##             expect_true(!is.null(reads.dt$anchor))
##             expect_true(!is.null(reads.dt$concord))
##             minus strand reads are reverse complemented
##             expect_true(all(reads.dt[strand == "-", seq != reading.frame]))
##             expect_true(all(reads.dt[strand == "+", seq == reading.frame]))
##         }
##     )
## })

## test_that(desc = "check that loose reads are put into the correct format", code = {
##     suppressWarnings(
##         expr = {
##             reads.dt = prep_loose_reads(li = le.dt, loose.reads.dt = loosereads.dt)
##             check that concord and anchor are added
##             expect_true(!is.null(reads.dt$anchor))
##             expect_true(!is.null(reads.dt$concord))
##             minus strand reads are reverse complemented
##             expect_true(all(reads.dt[strand == "-", seq != reading.frame]))
##             expect_true(all(reads.dt[strand == "+", seq == reading.frame]))
##         }
##     )
## })

## test_that(desc = "check that caller produces the expected call", code = {
##     suppressWarnings(
##         expr = {
##             call.res = call_loose_end(li = le.dt, ri = reads.dt,
##                                       concat.bwa = concat,
##                                       human.bwa = human,
##                                       mix.tn = TRUE,
##                                       verbose = FALSE)
##             expect_true(!is.null(call.res$call))
##             expect_true(!is.null(call.res$filtered.contigs))
##             expect_true(call.res$call$category == "type 1 loose end")
##             expect_true(call.res$filtered.contigs[, .N] > 0)
##         }
##     )
## })

test_that(desc = "check that caller with minimap produces expected call", code = {
    suppressWarnings(
        expr = {
            call.res = call_loose_end(li = le.dt, ri = reads.dt,
                                      concat.bwa = concat,
                                      human.bwa = human,
                                      concat.fn = concat.fasta.fn,
                                      pad = 5e3,
                                      max.reads = 2e4,
                                      outdir = "~/testing_tmp",
                                      minimap = TRUE,
                                      mix.tn = TRUE,
                                      verbose = FALSE)
            expect_true(!is.null(call.res$call))
            expect_true(!is.null(call.res$filtered.contigs))
            expect_true(call.res$call$category == "type 1 loose end")
            expect_true(call.res$filtered.contigs[, .N] > 0)
        }
    )
})


## test_that(desc = "test wrapper for caller", code = {
##     suppressWarnings(
##         expr = {
##             call.res = suppressWarnings(call_loose_end_wrapper(id = this.pair,
##                                                                le.dt = big.le.dt,
##                                                                reads.dt = big.reads.dt,
##                                                                concat.bwa = concat,
##                                                                human.bwa = human,
##                                                                pad = 1000,
##                                                                mix.tn = TRUE,
##                                                                minimap = FALSE,
##                                                                max.reads = 2e4,
##                                                                verbose = FALSE))
##             expect_true(!is.null(call.res$call))
##             expect_true(call.res$call[, .N] == 2)
##             expect_true(!is.null(call.res$filtered.contigs))
##             expect_true("type 1 loose end" %in% call.res$call$category)
##             expect_true(call.res$filtered.contigs[, .N] > 0)
##         }
##     )
## })

## test_that(desc = "test that caller wrapper works with minimap", code = {
##     suppressWarnings(
##         expr = {
##             call.res = suppressWarnings(call_loose_end_wrapper(id = this.pair,
##                                                                le.dt = big.le.dt,
##                                                                reads.dt = big.reads.dt,
##                                                                concat.bwa = concat,
##                                                                human.bwa = human,
##                                                                concat.fn = concat.fasta.fn
##                                                                pad = 5e3,
##                                                                mix.tn = TRUE,
##                                                                minimap = TRUE,
##                                                                outdir = "~/testing_tmp",
##                                                                max.reads = 2e4,
##                                                                verbose = FALSE))
##             expect_true(!is.null(call.res$call))
##             expect_true(call.res$call[, .N] == 2)
##             expect_true(!is.null(call.res$filtered.contigs))
##             expect_true("type 1 loose end" %in% call.res$call$category)
##             expect_true(call.res$filtered.contigs[, .N] > 0)
##         }
##     )
## })
