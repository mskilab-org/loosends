library(loosends)
library(testthat)

## reference directory
ref.dir = "~/git/loosends/inst/extdata/hg19_loosends"
loosereads.fn = system.file("tests", "testthat", "data", "loosereads_1", "loosereads.res.rds", package = "loosends")
reads.dt.fn = system.file("tests", "testthat", "data", "loosereads_1", "reads.dt.rds", package = "loosends")
le.dt.fn = system.file("tests", "testthat", "data", "loosereads_1", "le.dt.rds", package = "loosends")

## useful params
this.pair = "G32831.HCC1954"
this.le = parse.gr("20:60158837-")

## read input files
loosereads.dt = readRDS(loosereads.fn)
le.dt = readRDS(le.dt.fn)
reads.dt = readRDS(reads.dt.fn)

## load ref bwa
ref.obj = grab_ref_obj(ref.dir = ref.dir)

test_that(desc = "check that loose ends are put into correct format", code = {
    suppressWarnings(
        expr = {
            le.dt = prep_loose_ends(this.le, id = this.pair)
            expect_true(!is.null(le.dt$leix))
            expect_true(all(le.dt$sample == this.pair))
        }
    )
})

test_that(desc = "check that loose reads are put into the correct format", code = {
    suppressWarnings(
        expr = {
            reads.dt = prep_loose_reads(li = le.dt, loose.reads.dt = loosereads.dt)
            ## check that concord and anchor are added
            expect_true(!is.null(reads.dt$anchor))
            expect_true(!is.null(reads.dt$concord))
            ## minus strand reads are reverse complemented
            expect_true(all(reads.dt[strand == "-", seq != reading.frame]))
            expect_true(all(reads.dt[strand == "+", seq == reading.frame]))
        }
    )
})

test_that(desc = "check that loose reads are put into the correct format", code = {
    suppressWarnings(
        expr = {
            reads.dt = prep_loose_reads(li = le.dt, loose.reads.dt = loosereads.dt)
            ## check that concord and anchor are added
            expect_true(!is.null(reads.dt$anchor))
            expect_true(!is.null(reads.dt$concord))
            ## minus strand reads are reverse complemented
            expect_true(all(reads.dt[strand == "-", seq != reading.frame]))
            expect_true(all(reads.dt[strand == "+", seq == reading.frame]))
        }
    )
})

test_that(desc = "check that caller produces the expected call", code = {
    suppressWarnings(
        expr = {
            call.res = call_loose_end(li = le.dt, ri = reads.dt,
                                      ref_obj = ref.obj,
                                      mix.tn = TRUE,
                                      verbose = FALSE)
            expect_true(!is.null(call.res$call))
            expect_true(!is.null(call.res$filtered.contigs))
            expect_true(call.res$call$category == "type 1 loose end")
            expect_true(call.res$filtered.contigs[, .N] > 0)
        }
    )
})
