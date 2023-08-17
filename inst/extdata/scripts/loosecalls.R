{
    library(optparse)
    library(devtools)
    
    ## DO NOT FAIL SILENTLY
    ## options(error = function() {traceback(2); quit("no", 1)})

    if (!exists('opt')) {
        option_list = list(
            make_option(c("-i", "--id"), type = "character", help = "Sample / pair ID (required)"),
            make_option("--libdir", type = "character", help = "libdir", default=""),
            make_option("--outdir", type = "character", help = "outdir", default="./"),
            make_option("--ref_dir", type = "character", help = "ref dir", default="./"),
            make_option("--pad", type = "numeric", help = "assembly window (bp)", default=5e3),
            make_option("--max.qnames", type = "numeric", help = "max number of qnames", default=1e4),
            make_option("--mix", type = "logical", help = "mix t/n reads before assembly", default=FALSE),
            make_option("--overwrite", type = "logical", help = "overwrite existing analysis", default=TRUE),
            make_option("--minimap", type = "logical", help = "use minimap?", default=FALSE),
            make_option("--maxwin", type = "numeric", help = "max window around loose end to grab reads", default=5e3),
            make_option("--reads", type = "character", help = "data table of reads", default = "/dev/null"),
            make_option("--recall", type = "logical", help = "re-call?", default = TRUE),
            make_option("--cores", type = "numeric", help = "number of cores", default = 8),
            make_option("--loose", type = "character", help = "data table of loose ends", default="/dev/null"),
            make_option("--stash", type = "logical", help = "stash results and read from disk?", default=TRUE)
        )
        parseobj = OptionParser(option_list=option_list)
        opt = parse_args(parseobj)

        print(opt)
        
        print(.libPaths())
        options(error=function()traceback(2))

        ## keep record of run
        writeLines(paste(paste('--', names(opt), ' ', sapply(opt, function(x) paste(x, collapse = ',')), sep = '', collapse = ' '), sep = ''), paste(opt$outdir, 'cmd.args', sep = '/'))
        saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
    }

    message("Loading development version of loosends")
    devtools::load_all("~/git/loosends")
    library(gUtils)

    ## output filenames
    calls.dt.fn = file.path(opt$outdir, "calls.rds")
    recall.dt.fn = file.path(opt$outdir, "recalls.rds")
    contigs.dt.fn = file.path(opt$outdir, "contigs.rds")
    filtered.dt.fn = file.path(opt$outdir, "filtered.rds")


    ## read loose end file
    if (grepl(".rds", opt$loose, ignore.case = TRUE)) {
        loose.dt = readRDS(opt$loose)
        if (loose.dt[, .N]) {
            if (!("sample" %in% colnames(loose.dt))) {
                loose.dt[, sample := opt$id]
            }
            loose.dt = loose.dt[seqnames %in% c(as.character(1:22), "X", "Y"),]
        }
    } else {
        loose.dt = fread(opt$loose)[sample == opt$id,]
    }

    if (!nrow(loose.dt)) {
        message("No loose ends for this sample!")

        call.res = data.table()
        filtered.res = data.table()
        contig.res = data.table()
    } else {
        ## add required metadata
        message("Prepping loose ends")
        
        ## ALERT!!!!
        loose.dt = loose.dt[, sample := opt$id]

        ## add loose.end annotation if mising
        if (is.null(loose.dt$loose.end)) {
            loose.dt[, loose.end := paste0(seqnames, ":", start, strand)]
        }
        loose.dt = prep_loose_ends(li = loose.dt, id = opt$id)

        ## grab loose reads
        message("Grabbing reads")
        if (grepl(".rds", opt$reads, ignore.case = TRUE)) {
            loose.reads.dt = readRDS(opt$reads)
        } else {
            loose.reads.dt = fread(opt$reads)
        }

        ## convert reads to GRanges
        if (inherits(loose.reads.dt, 'data.table')) {
            loose.reads.gr = gUtils::dt2gr(loose.reads.dt[, .(seqnames, start, end)])
            loose.ends.gr = gUtils::dt2gr(loose.dt[, .(seqnames, start, end, sample, leix)])
        } else {
            stop("loose reads must be data.table")
        }

        ## load reference from disk
        message("Loading references from disk")

        ## only run if filtered contigs don't already exist
        if (TRUE) {
        ## if (opt$overwrite || (!file.exists(filtered.dt.fn))) {

            human.ref = paste0(opt$ref_dir, "/human_g1k_v37_decoy.fasta")
            concat.ref = paste0(opt$ref_dir, "/concatenated_references_deduped.fasta")

            human.bwa = RSeqLib::BWA(human.ref)
            concat.bwa = RSeqLib::BWA(concat.ref)

            ## if running minimap, pre-index the fasta to reduce some overhead
            if (opt$minimap) {
                message("Indexing minimap reference!")
                mmi = loosends::minimap_index(fasta = concat.ref, outdir = opt$outdir, verbose = TRUE)
            } else {
                message("Skipping indexing! Not using minimap.")
                mmi = concat.ref
            }

            ## divide up the reads by loose end ahead of time
            reads.per.le = lapply(1:loose.dt[, .N],
                                  function(ix, reads, loose, reads.gr) {
                                      message("Grabbing reads for loose end ", ix, " of ", loose[, .N])
                                      distances = GenomicRanges::distance(parse.gr(loose[ix, loose.end]),
                                                                          reads.gr,
                                                                          ignore.strand = TRUE)
                                      qnames = reads[which(distances < opt$maxwin),qname]
                                      n.qnames = length(unique(qnames))
                                      message("Number of unique qnames: ", n.qnames)
                                      if (n.qnames > opt$max.qnames) {
                                          message("Downsampling to ", opt$max.qnames, " qnames")
                                          qnames = sample(qnames, size = opt$max.qnames, replace = FALSE)
                                      }
                                      message("Prepping reads...")
                                      this.loose.reads = reads[qname %in% qnames,]
                                      ri = prep_loose_reads(loose[ix,], this.loose.reads)
                                      return(ri)
                                  },
                                  loose.reads.dt,
                                  loose.dt,
                                  loose.reads.gr)

            ## split loose ends into lsit
            le.list = lapply(1:loose.dt[, .N],
                             function(ix, loose) {
                                 return(loose[ix,])
                             },
                             loose.dt)

            ## remove stuff to save memory
            ## cuz reads can be big
            rm(loose.reads.dt)
            rm(loose.reads.gr)
            gc()

            ## get the size of reads
            ## if it is more than 0.5 gb don't fork so many times??
            if (file.info(opt$reads)$size > 5e8) {
                message("large input file detected! using fewer cores to prevent memory issues from forking.")
                opt$cores = 2
            }

            if (opt$minimap) {
                message("using minimap... going single core mode :(")
                opt$cores = 1
            }

            if (opt$stash) {
                if (!dir.exists(file.path(opt$outdir, "tmp_loosends"))) {
                    message("Creating directory for stashing outputs")
                    dir.create(file.path(opt$outdir, "tmp_loosends"), recursive = TRUE)
                }
            }

            message("Starting caller")
            browser()
            res = mcmapply(FUN = function(li, ri, human, concat) {
                message("Starting: ", paste0(li[, seqnames], ":", li[, start], li[, strand]))
                ## create file to shash results
                stash.fn = file.path(opt$outdir, "tmp_loosends", paste0(li[, leix], ".rds"))
                ## create directory to stash minimap stuff
                stash.dir = file.path(opt$outdir, "tmp_loosends", li[, leix])
                dir.create(path = stash.dir, recursive = TRUE)
                if (opt$stash && file.exists(stash.fn) && (file.info(stash.fn)$size > 0)) {
                    message("Reading results from stash")
                    call.res = readRDS(stash.fn)
                    message("Annotation: ", call.res$label.dt$annotation)
                    message("Somatic: ", call.res$label.dt$somatic)
                    return(call.res)
                } else {
                    call.res = loosends:::call_loose_end2(li,
                                                          ri,
                                                          id = opt$id,
                                                          concat.bwa = concat,
                                                          human.bwa = human,
                                                          window = opt$pad,
                                                          use.minimap = opt$minimap,
                                                          outdir = stash.dir,
                                                          verbose = (opt$cores == 1))
                    if (opt$stash) {
                        saveRDS(call.res, stash.fn)
                    }
                }
                message("Annotation: ", call.res$label.dt$annotation)
                message("Somatic: ", call.res$label.dt$somatic)
                return(call.res)
            },
            le.list,
            reads.per.le,
            MoreArgs = list(human = human.bwa, concat = concat.bwa),
            mc.cores = opt$cores,
            mc.preschedule = FALSE,
            SIMPLIFY = FALSE)

            message("Concatenating calls")
            call.res = lapply(res, function(x) {return(x$label.dt)}) %>%
                rbindlist(fill = TRUE)
            contig.res = lapply(res, function(x) {return(x$all.contigs)}) %>%
                rbindlist(fill = TRUE)
            filtered.res = lapply(res, function(x) {return(x$keep.tigs.support)}) %>%
                rbindlist(fill = TRUE)

            call.res[, sample := opt$id]

            if (filtered.res[, .N]) {
                filtered.res[, sample := opt$id]
            }
            message("Saving output")
            saveRDS(call.res, calls.dt.fn)
            saveRDS(contig.res, contigs.dt.fn)
            saveRDS(filtered.res, filtered.dt.fn)
        } else {
            message("Filtered contigs already exist!!")
            filtered.res = readRDS(filtered.dt.fn)
        }
    }

    message("Done")

    quit("no")
}
