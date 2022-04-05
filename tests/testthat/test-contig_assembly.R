library(loosends)
library(gUtils)

dup.calns.fn = system.file("extdata", "tests", "new_caller_2", "dup.calns.rds", package = "loosends")
del.calns.fn = system.file("extdata", "tests", "new_caller_2", "del.calns.rds", package = "loosends")
inv.calns.fn = system.file("extdata", "tests", "new_caller_2", "inv.calns.rds", package = "loosends")
fbi.calns.fn = system.file("extdata", "tests", "new_caller_2", "fbi.calns.rds", package = "loosends")
tra.calns.fn = system.file("extdata", "tests", "new_caller_2", "tra.calns.rds", package = "loosends")
cpx.calns.fn = system.file("extdata", "tests", "new_caller_2", "cpx.calns.rds", package = "loosends")

dup.calns = readRDS(dup.calns.fn)
del.calns = readRDS(del.calns.fn)
inv.calns = readRDS(inv.calns.fn)
fbi.calns = readRDS(fbi.calns.fn)
tra.calns = readRDS(tra.calns.fn)
cpx.calns = readRDS(cpx.calns.fn)

## SEEDS
dup.grl = parse.grl("16:30070288-30070288+,16:30079541-30079541-")
del.grl = parse.grl("16:30070288-30070288-,16:30079541-30079541+")
inv.grl = parse.grl("4:75231865-75231865-,4:75233388-75233388-")
fbi.grl = parse.grl("4:75231865-75231865-,4:75231865-75231865-")
tra.grl = parse.grl("10:110480091-110480091-,4:167114710-167114710+")
cpx.grl = parse.grl("10:110480091-110480091-,4:167114710-167114710+,16:30079541-30079541-")

## get the seed region, just the first range
dup.seed = gr.flipstrand(GenomicRanges::resize(dup.grl[[1]][1], width = 250, fix = "start"))
del.seed = gr.flipstrand(GenomicRanges::resize(del.grl[[1]][1], width = 250, fix = "start"))
inv.seed = gr.flipstrand(GenomicRanges::resize(inv.grl[[1]][1], width = 250, fix = "start"))
fbi.seed = gr.flipstrand(GenomicRanges::resize(fbi.grl[[1]][1], width = 250, fix = "start"))
tra.seed = gr.flipstrand(GenomicRanges::resize(tra.grl[[1]][1], width = 250, fix = "start"))
cpx.seed = gr.flipstrand(GenomicRanges::resize(cpx.grl[[1]][1], width = 250, fix = "start"))


## REFERNCE GRANGES
unassembled.gr = readRDS(system.file("extdata", "assembly.gaps.ranges.rds"))
low.mappability.gr = readRDS(system.file("extdata", "wide.multimapping.ranges.rds"))

test_that(desc = "test contig QC for DUP", code = {
    suppressWarnings(
        expr = {
            dup.calns = add_contig_ctypes(dup.calns)
            dup.qc = qc_contigs(calns = as.data.table(dup.calns),
                                seed.gr = dup.seed,
                                unassembled.gr = unassembled.gr,
                                low.mappability.gr = low.mappability.gr)

            ## check classification is correct
            expect_true(all(dup.qc[, junction]))
            expect_true(all(dup.qc[, !is.na(jstring)]))

            ## check correct jstring
            gr = unlist(parse.grl(dup.qc[, jstring][1]))
            expect_true(all(gr %^% (unlist(dup.grl) + 250)))

            ## check orientation is correct
            expect_true(strand(gr[1]) == "+")
            expect_true(strand(gr[2]) == "-")
        })
})

test_that(desc = "test contig QC for DEL", code = {
    suppressWarnings(
        expr = {
            del.calns = add_contig_ctypes(del.calns)
            del.qc = qc_contigs(calns = del.calns,
                                seed.gr = del.seed,
                                unassembled.gr = unassembled.gr,
                                low.mappability.gr = low.mappability.gr)

            ## check classification is correct
            expect_true(all(del.qc[, junction]))
            expect_true(all(del.qc[, !is.na(jstring)]))

            ## check correct jstring
            gr = unlist(parse.grl(del.qc[, jstring][1]))
            expect_true(all(gr %^% (unlist(del.grl) + 250)))

            ## check orientation is correct
            expect_true(strand(gr[1]) == "-")
            expect_true(strand(gr[2]) == "+")
        })
})

test_that(desc = "test contig QC for TRA", code = {
    suppressWarnings(
        expr = {
            tra.calns = add_contig_ctypes(tra.calns)
            tra.qc = qc_contigs(calns = as.data.table(tra.calns),
                                seed.gr = tra.seed,
                                unassembled.gr = unassembled.gr,
                                low.mappability.gr = low.mappability.gr)

            ## check classification is correct
            expect_true(all(tra.qc[, junction]))
            expect_true(all(tra.qc[, !is.na(jstring)]))

            ## check correct jstring
            gr = unlist(parse.grl(tra.qc[, jstring][1]))
            expect_true(all(gr %^% (unlist(tra.grl) + 250)))

            ## check orientation is correct
            expect_true(strand(gr[1]) == "-")
            expect_true(strand(gr[2]) == "+")
        })
})

test_that(desc = "test contig QC for INV", code = {
    suppressWarnings(
        expr = {
            inv.calns = add_contig_ctypes(inv.calns)
            inv.qc = qc_contigs(calns = as.data.table(inv.calns),
                                seed.gr = inv.seed,
                                unassembled.gr = unassembled.gr,
                                low.mappability.gr = low.mappability.gr)

            ## check classification is correct
            expect_true(all(inv.qc[, junction]))
            expect_true(all(inv.qc[, !is.na(jstring)]))

            ## check correct jstring
            gr = unlist(parse.grl(inv.qc[, jstring][1]))
            expect_true(all(gr %^% (unlist(inv.grl) + 250)))

            ## check orientation is correct
            expect_true(strand(gr[1]) == "-")
            expect_true(strand(gr[2]) == "-")
        })
})

test_that(desc = "test contig QC for FBI", code = {
    suppressWarnings(
        expr = {
            fbi.calns = add_contig_ctypes(fbi.calns)
            fbi.qc = qc_contigs(calns = as.data.table(fbi.calns),
                                seed.gr = fbi.seed,
                                unassembled.gr = unassembled.gr,
                                low.mappability.gr = low.mappability.gr)

            ## check classification is correct
            expect_true(all(fbi.qc[, junction]))
            expect_true(all(fbi.qc[, fbi]))
            expect_true(all(fbi.qc[, !is.na(jstring)]))

            ## check correct jstring
            gr = unlist(parse.grl(fbi.qc[, jstring][1]))
            expect_true(all(gr %^% (unlist(fbi.grl) + 250)))

            ## check orientation is correct
            expect_true(strand(gr[1]) == "-")
            expect_true(strand(gr[2]) == "-")
        })
})

test_that(desc = "test contig QC for CPX", code = {
    suppressWarnings(
        expr = {
            cpx.calns = add_contig_ctypes(cpx.calns)
            cpx.qc = qc_contigs(calns = as.data.table(cpx.calns),
                                seed.gr = cpx.seed,
                                unassembled.gr = unassembled.gr,
                                low.mappability.gr = low.mappability.gr)

            ## check classification is correct
            expect_true(all(fbi.qc[, junction]))
            expect_true(all(fbi.qc[, fbi]))
            expect_true(all(fbi.qc[, !is.na(jstring)]))

            ## check correct jstring
            gr = unlist(parse.grl(fbi.qc[, jstring][1]))
            expect_true(all(gr %^% (unlist(fbi.grl) + 250)))

            ## check orientation is correct
            expect_true(strand(gr[1]) == "-")
            expect_true(strand(gr[2]) == "-")
        })
})
