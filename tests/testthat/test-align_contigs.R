library(loosends)
library(gUtils)

concat.fn = system.file("extdata", "hg19_loosends", "concatenated_references_deduped.fasta", package = "loosends")
concat = BWA(fasta = concat.fn)

dup.tigs.fn = system.file("extdata", "tests", "new_caller_1", "dup.tigs.rds", package = "loosends")
inv.tigs.fn = system.file("extdata", "tests", "new_caller_1", "inv.tigs.rds", package = "loosends")
tra.tigs.fn = system.file("extdata", "tests", "new_caller_1", "tra.tigs.rds", package = "loosends")

dup.tigs = readRDS(dup.tigs.fn)
inv.tigs = readRDS(inv.tigs.fn)
tra.tigs = readRDS(tra.tigs.fn)

## hard code the expected outputs
dup.grl = parse.grl("16:30079541-30079541-,16:30070288-30070288+")
inv.grl = parse.grl("4:75231865-75231865-,4:75233388-75233388-")
tra.grl = parse.grl("10:110480091-110480091-,4:167114710-167114710+")

## get breakpoints
dup.bp = stack(dup.grl)
inv.bp = stack(inv.grl)
tra.bp = stack(tra.grl)

test_that(desc = "check that alignment finds the correct breakends", code = {
    suppressWarnings(
        expr = {
            dup.aln = align_contigs(tigs = dup.tigs, ref = concat)
            expect_true(all((dup.bp + 1e3) %^% dt2gr(dup.aln)))
            inv.aln = align_contigs(tigs = inv.tigs, ref = concat)
            expect_true(all((inv.bp + 1e3) %^% dt2gr(inv.aln)))
            tra.aln = align_contigs(tigs = tra.tigs, ref = concat)
            expect_true(all((tra.bp + 1e3) %^% dt2gr(tra.aln)))
        })
})

dup.aln.fn = system.file("extdata", "tests", "new_caller_1", "dup.aln.rds", package = "loosends")
inv.aln.fn = system.file("extdata", "tests", "new_caller_1", "inv.aln.rds", package = "loosends")
tra.aln.fn = system.file("extdata", "tests", "new_caller_1", "tra.aln.rds", package = "loosends")

inv.bp.gr.fn = system.file("extdata", "tests", "new_caller_1", "inv.bp.gr.rds", package = "loosends")
inv.bp.gr = readRDS(inv.bp.gr.fn)
dup.bp.gr.fn = system.file("extdata", "tests", "new_caller_1", "dup.bp.gr.rds", package = "loosends")
dup.bp.gr = readRDS(dup.bp.gr.fn)
tra.bp.gr.fn = system.file("extdata", "tests", "new_caller_1", "tra.bp.gr.rds", package = "loosends")
tra.bp.gr = readRDS(tra.bp.gr.fn)

dup.aln = readRDS(dup.aln.fn)
inv.aln = readRDS(inv.aln.fn)
tra.aln = readRDS(tra.aln.fn)

test_that(desc = "check that contig QC retains informative contigs", code = {
    suppressWarnings(
        expr = {
            dup.aln.qc = qc_contigs(dup.aln, gr.flipstrand(dup.bp.gr) + 1e3)
            expect_true(dup.aln.qc[(keep), .N] > 0)
            expect_true(dup.aln.qc[(keep) & (outside.seed), .N] > 0)
            
            inv.aln.qc = qc_contigs(inv.aln, gr.flipstrand(inv.bp.gr) + 1e3)
            expect_true(inv.aln.qc[(keep), .N] > 0)
            expect_true(inv.aln.qc[(keep) & (fbi), .N] > 0) ## check for FBI specifically

            tra.aln.qc = qc_contigs(tra.aln, gr.flipstrand(tra.bp.gr) + 1e3)
            expect_true(tra.aln.qc[(keep), .N] > 0)
            expect_true(tra.aln.qc[(keep) & (outside.seed), .N] > 0)
            expect_true(any(dt2gr(tra.aln.qc[(keep), .(seqnames, start, end, strand)]) %^% tra.bp))
        })
})

