{
    library(optparse)
    library(devtools)
    library(gUtils)

    ## do not fail silently!
    options(error = function() {traceback(2); quit("no", 1)})

    if (!exists('opt')) {
        option_list = list(
            make_option(c("-i", "--id"), type = "character", help = "Sample / pair ID (required)"),
            make_option("--tbam", type = "character", help = "Path to tumor bam (required)"),
            make_option("--nbam", type = "character", help = "Path to normal bam (required)"),
            make_option("--ref", type = "character", help = "path to fasta for realignment (required)"),
            make_option("--loose", type = "character",
                        help = "path to loose ends file - either gGraph/JaBbA output or text file (required)"),
            make_option("--overwrite", type = "logical",
                        help = "overwrite existing analyses in output directory?",
                        default=TRUE),
            make_option("--pad", type = "numeric",
                        help = "loose end window padding",
                        default = 5e3),
            make_option("--anchorlift_window", type = "numeric",
                        help = "loose read anchorlift window",
                        default = 1e4),
            make_option("--bowtie", type = "logical", help = "use bowtie?", default=FALSE),
            make_option("--bowtie_dir", type = "character", help = "path to bowtie index directory",
                        default="~/git/loosends/inst/extdata/hg19_loosends"),
            make_option("--bowtie_ref", type = "character", help = "basename of bowtie index files",
                        default="human_g1k_v37_decoy"),
            make_option("--libdir", type = "character", help = "libdir", default=""),
            make_option("--outdir", type = "character", help = "outdir", default="./")
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
    library(gGnome)

    ## output filenames
    loose.ends.fn = paste0(opt$outdir, "/", "loose.ends.rds")
    all.reads.fn = paste0(opt$outdir, "/", "all.reads.rds")
    anchor.fn = paste0(opt$outdir, "/", "anchor.reads.rds")
    

    if (opt$overwrite | !file.exists(loose.ends.fn)) {
        ## read loose end file
        if (grepl(".rds", opt$loose, ignore.case = TRUE)) {

            loose.dt = readRDS(opt$loose)

            ## check for gGraph input
            if (inherits(loose.dt, 'list')) {
                if (all(c('segstats', 'adj', 'purity', 'ploidy') %in% names(loose.dt))) {
                    message("Detected JaBbA output!")
                    jab = gG(jabba = loose.dt)
                    loose.dt = as.data.table((jab$loose %Q% (terminal == FALSE)))
                } else {
                    stop("Invalid file supplied")
                }
            } else if (inherits(loose.dt, 'gGraph')) {
                message("Detected gGraph input!")
                loose.dt = as.data.table((loose.dt$loose %Q% (terminal == FALSE)))
            } else if (inherits(loose.dt, 'data.table')) {
                message("Detected data.table!")
            } else if (inherits(loose.dt, 'data.frame')) {
                message("Detected data.frame")
                loose.dt = as.data.table(loose.dt)
            } else if (inherits(loose.dt, 'GRanges')) {
                loose.dt = as.data.table(loose.dt)
            } else {
                stop("Invalid file supplied for loose ends")
            }
            
        } else {
            loose.dt = fread(opt$loose)[sample == opt$id,]
        }

        ## check for nonempty table
        if (loose.dt[, .N]) {
            if (!("sample" %in% colnames(loose.dt))) {
                loose.dt[, sample := opt$id]
            }
            loose.dt = loose.dt[(seqnames %in% c(as.character(1:22), "X", "Y")) |
                                (seqnames %in% paste0('chr', c(as.character(1:22), "X", "Y"))), ]
        }

        message("Saving loose ends for this sample...")
        saveRDS(loose.dt, loose.ends.fn)
    } else {
        loose.dt = readRDS(loose.ends.fn)
    }


    ## if there are no loose ends, save empty data tables and exit
    if (!nrow(loose.dt)) {
        message("Sample has no loose ends!")
        saveRDS(data.table(), all.reads.fn)
        saveRDS(data.table(), anchor.fn)
    } else {

        if (opt$overwrite | !file.exists(all.reads.fn) | file.info(all.reads.fn)$size == 0) {
            message("Grabbing loose reads")

            if (opt$bowtie) {
                reads.dt = loosereads_wrapper(dt2gr(loose.dt),
                                              tbam = opt$tbam,
                                              nbam = opt$nbam,
                                              ref = opt$ref,
                                              id = opt$id,
                                              outdir = opt$outdir,
                                              bowtie = TRUE,
                                              bowtie.ref = opt$bowtie_ref,
                                              bowtie.dir = opt$bowtie_dir,
                                              pad = opt$pad,
                                              verbose = TRUE)
            } else {
                reads.dt = loosereads_wrapper(dt2gr(loose.dt),
                                              tbam = opt$tbam,
                                              nbam = opt$nbam,
                                              ref = opt$ref,
                                              id = opt$id,
                                              outdir = opt$outdir,
                                              pad = opt$pad,
                                              verbose = TRUE)
            }
            saveRDS(reads.dt, all.reads.fn)
        } else {
            message("Reading loose reads from file")
            reads.dt = readRDS(all.reads.fn)
        }

        if (opt$overwrite | !file.exists(anchor.fn) | file.info(anchor.fn)$size == 0) {
            if (!nrow(reads.dt)) {
                message("No loose reads!")
                saveRDS(data.table(), anchor.fn)
            } else {

                message("Annotating loose pair seeds")
                anchor.tmp.dt = reads.dt[(high.mate) & (loose.pair),]
                if (!nrow(anchor.tmp.dt)) {
                    anchor.dt = data.table()
                } else {
                    anchor.gr = dt2gr(anchor.tmp.dt[, .(seqnames, start, end, strand)],
                                      seqlengths = hg_seqlengths())
                    loose.gr = dt2gr(loose.dt[, .(seqnames, start, end, strand)],
                                     seqlengths = hg_seqlengths())

                    ## actually use gUtils anchorlift
                    anchorlift.gr = gUtils::anchorlift(query = anchor.gr,
                                                       subject = loose.gr,
                                                       window = opt$anchorlift_window)

                    anchor.dt = as.data.table(anchorlift.gr)

                    
                    ## add annotations
                    anchor.dt[, sample := anchor.tmp.dt$sample[query.id]]
                    anchor.dt[, track := anchor.tmp.dt$track[query.id]]
                    anchor.dt[, tumor := ifelse(track %like% "sample", "tumor", "normal")]
                    anchor.dt[, query.strand := anchor.tmp.dt$strand[query.id]]
                    anchor.dt[, subject.strand := loose.dt$strand[subject.id]]

                    ## flip the strands because the loose ends strands are in "junction orientation"
                    anchor.dt[, forward := query.strand != subject.strand]

                    ## add whether the loose end is annotated as 'germline' or 'somatic'
                    anchor.dt[, somatic := NA]
                    if (!is.null(loose.dt$keep)) {
                        anchor.dt[, somatic := ifelse(loose.dt$keep[subject.id], "somatic", "germline")]
                    }
                    
                }

                saveRDS(anchor.dt, anchor.fn)
            }
        } else {
            message("anchorlift file already exists!")
        }
        
        message("Done!")
    }
    quit("no", status = 0)
}
