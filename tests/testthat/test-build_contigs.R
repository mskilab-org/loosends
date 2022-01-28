library(loosends)

## check that contigs can be assembled from pre-saved reads
dup.seed.rds.fn = system.file("extdata", "tests", "new_caller_1", "dup.seed.rds", package = "loosends")
tra.seed.rds.fn = system.file("extdata", "tests", "new_caller_1", "tra.seed.rds", package = "loosends")
inv.seed.rds.fn = system.file("extdata", "tests", "new_caller_1", "inv.seed.rds", package = "loosends")

dup.seed.rds = readRDS(dup.seed.rds.fn)
tra.seed.rds = readRDS(tra.seed.rds.fn)
inv.seed.rds = readRDS(inv.seed.rds.fn)

test_that(desc = "check that assembly will find contigs at junction breakends", code = {
    suppressWarnings(
        expr = {
            tra.tigs = build_contigs(tra.seed.rds)
            expect_true(length(tra.tigs) > 0)
            dup.tigs = build_contigs(dup.seed.rds)
            expect_true(length(dup.tigs) > 0)
            inv.tigs = build_contigs(inv.seed.rds)
            expect_true(length(inv.tigs) > 0)
        })
})
