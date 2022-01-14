#' @import GenomicRanges
#' @import data.table
#' @import Matrix
#' @import igraph
#' @import Rsamtools
#' @import reshape2
#' @import rtracklayer
#' @import gUtils
#' @import bamUtils
#' @import Biostrings
#' @import RSeqLib
#' @import gChain
#' @import gTrack
#' @import gGnome
#' @import BSgenome.Hsapiens.UCSC.hg19

#' @name loosereads_wrapper
#' @title loosereads_wrapper
#'
#' @description
#'
#' wrapper for obtaining reads + mates around loose ends and realigning
#'
#' @param ranges (GRanges) loose end GRanges
#' @param tbam (character) path to normal bam
#' @param nbam (character) optional, path to normal bam
#' @param ref (character) path to reference .fasta
#' @param bowtie (logical) use bowtie for realignment?
#' @param bowtie.ref (character) basename of bowtie references
#' @param bowtie.dir (character) dirname of bowtie
#' @param id (character) sample name
#' @param outdir (character)
#' @param pad (numeric)
#' @param verbose (logical)
#'
#' @return data.table of reads + mates
loosereads_wrapper = function(ranges = GRanges(),
                              tbam = "/dev/null",
                              nbam = "/dev/null",
                              ref = "/dev/null",
                              bowtie = FALSE,
                              bowtie.ref = "/dev/null",
                              bowtie.dir = "/dev/null",
                              id = "",
                              outdir = "./",
                              pad = 5000,
                              verbose = FALSE) {

    tparams = grab_looseread_params(gr = ranges, bam = tbam, pad = pad, verbose = verbose)
    tsub = grab_loosereads(bam = tbam,
                           ranges = tparams$ranges,
                           qnames = tparams$qnames,
                           outdir = paste0(outdir, "/tumor"),
                           verbose = verbose)

    if (bowtie) {
        taln = realign_loosereads(bam = tsub,
                                  ref = bowtie.ref,
                                  bowtie = TRUE,
                                  bowtie.dir = bowtie.dir,
                                  outdir = paste0(outdir, "/tumor"),
                                  verbose = verbose)
    } else {
        taln = realign_loosereads(bam = tsub, ref = ref,
                                  outdir = paste0(outdir, "/tumor"),
                                  verbose = verbose)
    }
    
    if (check_file(nbam)) {
        nparams = grab_looseread_params(gr = ranges, bam = nbam, pad = pad, verbose = verbose)
        nsub = grab_loosereads(bam = nbam,
                               ranges = nparams$ranges,
                               qnames = nparams$qnames,
                               outdir = paste0(outdir, "/normal"),
                               verbose = verbose)

        if (bowtie) {
            naln = realign_loosereads(bam = nsub,
                                      ref = bowtie.ref,
                                      bowtie = TRUE,
                                      bowtie.dir = bowtie.dir,
                                      outdir = paste0(outdir, "/normal"),
                                      verbose = verbose)
        } else {
            naln = realign_loosereads(bam = nsub, ref = ref,
                                      outdir = paste0(outdir, "/normal"),
                                      verbose = verbose)
        }
        
    } else {
        nsub = "/dev/null"
        naln = "/dev/null"
    }

    res = loose.reads2(tbam = tsub, taln = taln, nbam = nsub, naln = naln, id = id, 
                       filter = FALSE, verbose = verbose)
    return(res)
        
}

#' @name check_file
#' @title check_file
#'
#' @description
#'
#' Check if file supplied is nonempty
#'
#' @param fn
#'
#' @return logical, TRUE if file exists and is nonempty
check_file = function(fn = NULL) {
    if (is.null(fn)) {
        return(FALSE)
    }

    if (is.na(fn)) {
        return(FALSE)
    }

    if (!file.exists(fn)) {
        return(FALSE)
    }

    if (!file.info(fn)$size) {
        return(FALSE)
    }

    return(TRUE)
}

#' @name has_chr
#' @title has_chr
#'
#' @description
#'
#' check if bam file has chr prefix
#'
#' @param bam
#'
#' @return logical, TRUE if bam seqnames start with chr
has_chr = function(bam) {

    if (!check_file(bam)) {
        stop("Invalid file supplied: ", bam)
    }

    bf = Rsamtools::BamFile(file = bam)
    res = any(grepl("chr", seqnames(seqinfo(bf))))
    return(res)
}

#' @name grab_looseread_params
#' @title grab_looseread_params
#'
#' @description
#'
#' Grab ranges and qnames of reads near loose ends + their mates
#' (Including ranges for split reads)
#'
#' @param gr (GRanges) stranded GRanges representing loose ends
#' @param bam (character) path to bam file
#' @param pad (numeric) pad around loose ends
#' @param mate_pad (numeric) pad around mate windows
#' @param verbose (logical)
#'
#' @return list with names:
#' - $ranges (GRanges of windows)
#' - $qnames (character vector of loose read qnames)
grab_looseread_params = function(gr = GRanges(),
                                 bam = character(),
                                 pad = 5000,
                                 mate_pad = 150,
                                 verbose = FALSE) {

    empty.res = list(ranges = GRanges(), qnames = character())

    if (!inherits(gr, "GRanges")) {
        stop("gr must be GRanges")
    }

    if (!length(gr)) {
        return(empty.res)
    }

    loose.gr = gr + pad

    if (has_chr(bam)) {
        if (verbose) {
            message("Adding chromosome prefix before reading bam")
        }
        loose.gr = gr.chr(loose.gr)
    }

    if (verbose) {
        message("Reading BAM: ", bam)
    }

    all.reads.grl = bamUtils::read.bam(bam = bam, intervals = loose.gr, pairs.grl = TRUE,
                                       isDuplicate = NA, isPaired = TRUE, tag = "SA")
    reads = as.data.table(unlist(all.reads.grl))

    if (!reads[, .N]) {
        if (verbose) {
            message("No reads in specified region!")
        }
        return(empty.res)
    }

    if (verbose) {
        message("Number of reads: ", reads[, .N])
        message("Checking windows of split alignments")
    }

    ## check for split reads and grab their locations
    splits = reads[!is.na(SA),]
    if (splits[, .N]) {
        splits$SA = as.character(splits$SA)
        splwin = dunlist(strsplit(splits$SA, ";"))
        spl = unlist(lapply(strsplit(splwin$V1, ","), function(w) paste(w[1], w[2], sep=":")))
        spl = GRanges(spl)
        spl.wins = trim(gUtils::gr.reduce(spl + mate_pad) - mate_pad)
    } else {
        spl.wins = GRanges()
    }

    if (verbose) {
        message("Number of split windows: ", length(spl.wins))
        message("Getting mate windows")
    }
    
    ## then grab all mate windows
    mw = reads[!is.na(mrnm) & !is.na(seq)]
    if (mw[, .N]) {
        mw[, ":="(seqnames = mrnm,
                  start = ifelse(is.na(mpos), start, mpos),
                  end = ifelse(is.na(mpos), end, mpos))]
        mw = dt2gr(mw[,!"strand"], seqlengths=seqlengths(loose.gr))
        mw[width(mw) == 0] = mw[width(mw) == 0] + 1
        mate.wins = gUtils::gr.reduce(mw+150)-150
    } else {
        mate.wins = GRanges(seqlengths = seqlengths(loose.gr))
    }

    if (verbose) {
        message("Number of mate windows: ", length(mate.wins))
    }

    windows = trim(reduce(grbind(loose.gr, mate.wins + mate_pad, spl.wins + mate_pad)))
    qnames = unique(reads[, qname])

    return(list(ranges = windows, qnames = qnames))
}

#' @name grab_loosereads
#' @title grab_loosereads
#'
#' @description
#'
#' get the reads + mates associated with loose ends
#'
#' @param bam (character) path to BAM file
#' @param ranges (GRanges) 
#' @param qnames (character)
#' @param outdir (character) path to dump results
#' @param overwrite (logical) overwrite existing files? default TRUE
#' @param verbose (logical) default FALSE
#'
#' @return path to indexed BAM file with all required reads
grab_loosereads = function(bam = NA_character_, ranges = GRanges(), qnames = character(),
                           outdir = "./", overwrite = TRUE, verbose = TRUE) {

    if (!dir.exists(outdir)) {
        dir.create(outdir, recursive = TRUE)
    }

    ## write windows and qnames
    regions.fn = file.path(outdir, "windows.bed")
    qnames.fn = file.path(outdir, "qnames.txt")

    if (verbose) {
        message("Saving regions bed file: ", regions.fn)
    }
    rtracklayer::export.bed(ranges, regions.fn)

    if (verbose) {
        message("Saving qnames: ", qnames.fn)
    }
    writeLines(qnames, con = qnames.fn, sep = "\n")

    noheader.fn = file.path(outdir, "noheader.sam")
    ## grab reads corresponding with desired qnames
    cmd = paste("samtools view -M -L", regions.fn, bam, "|",
                "LC_ALL=C grep -w -F -f", qnames.fn, ">",
                noheader.fn)

    if (verbose) {
        message("Grabbing loose reads...")
        message(cmd)
    }

    sys.res = system(cmd)
    if (sys.res) {
        file.remove(noheader.fn)
        stop("Error!")
    }

    header.fn = file.path(outdir, "header.sam")
    final.fn = file.path(outdir, "loosereads.bam")

    cmd = paste("samtools view -H ", bam, ">", header.fn)

    if (verbose) {
        message("Reheadering!")
        message(cmd)
    }


    sys.res = system(cmd)
    if (sys.res) {
        file.remove(header.fn)
        stop("Error!")
    }

    cmd = paste("cat", header.fn, noheader.fn, "|",
                "samtools view -S -b >", final.fn)

    if (verbose) {
        message(cmd)
    }
    
    sys.res = system(cmd)
    if (sys.res) {
        file.remove(final.fn)
        stop("Error!")
    }

    cmd = paste("samtools index", final.fn)

    if (verbose) {
        message(cmd)
    }
    
    sys.res = system(cmd)
    if (sys.res) {
        stop("Error!")
    }

    return(final.fn)
}

#' @name realign_loosereads
#' @title realign_loosereads
#'
#' @description
#'
#' Realign reads to selected reference
#'
#' @param bam (character) path to bam file
#' @param ref (character) path to BWA indexed reference. If running bowtie then this should be basename of the indexed reference file
#' @param bowtie (logical) use Bowtie2 aligner? default FALSE
#' @param bowtie.dir (logical) directory containing bowtie files
#' @param outdir (character) output directory to dump results
#' @param overwrite (logical) overwrite existing analysis
#' @param verbose (logical)
#'
#' @return path to output bam file with realigned reads
realign_loosereads = function(bam,
                              ref = system.file('extdata', 'hg19_looseends', 'human_g1k_v37_decoy.fasta', package='loosends'),
                              bowtie = FALSE,
                              bowtie.dir = system.file("extdata", "hg19_loosends", package = "loosends"),
                              outdir = "./",
                              overwrite = FALSE,
                              verbose = FALSE) {

    if (!dir.exists(outdir)) {
        dir.create(outdir, recursive = TRUE)
    }

    ## create FASTQ
    sorted.fn = file.path(outdir, "name.sorted.bam")
    fastq.fn = file.path(outdir, "loosereads.fq")

    cmd = paste("samtools sort -n ", bam, ">", sorted.fn)
    if (verbose) {
        message("Sorting reads by name")
        message(cmd)
    }
    sys.res = system(cmd)
    if (sys.res) {
        file.remove(sorted.fn)
        stop("Error!")
    }

    cmd = paste("samtools fastq -N ", sorted.fn, ">", fastq.fn)
    if (verbose) {
        message("Making fastq")
        message(cmd)
    }
    sys.res = system(cmd)
    if (sys.res) {
        file.remove(fastq.fn)
        stop("Error!")
    }

    aligned.fn = file.path(outdir, "aln.sam")
    final.fn = file.path(outdir, "aln.bam")

    if (bowtie) {

        ## use BOWTIE!!!
        if (verbose) {message("Using Bowtie2 for single end realignment")}

        ## change to outdir
        
        currdir  = normalizePath(getwd())
        fastq.fn = normalizePath(fastq.fn)
        aligned.fn = file.path(currdir, "aln.sam")

        if (verbose) {
            message("Current directory: ", currdir)
            message("Changing to bowtie directory for alignment: ", bowtie.dir)
        }
        setwd(bowtie.dir)

        ## cmd = paste("module load bowtie2; bowtie2 --very-sensitive-local -x", ref,
        ## bowtie2 should be runnable from command line but if not may need to load the module
        cmd = paste("bowtie2 --very-sensitive-local -x", ref,
                    "-U", fastq.fn,
                    "-S", aligned.fn)
        if (verbose) message("RUNNING: ", cmd)
        sys.res = system(cmd)
        if (sys.res) {
            file.remove(aligned.fn)
            stop("Error!")
        }

        if (verbose) {message("Changing back to output dir: ", currdir)}
        setwd(currdir)
    } else {

        ## perform single end realignment
        if (verbose) { message("Using BWA for single end realignment") }
        
        cmd = paste("bwa mem", ref, fastq.fn, ">", aligned.fn)
        if (verbose) {
            message("Realigning reads")
            message(cmd)
        }
        sys.res = system(cmd)
        if (sys.res) {
            file.remove(aligned.fn)
            stop("Error!")
        }
    }

    cmd = paste("samtools view -Sb", aligned.fn, "|",
                "samtools sort >", final.fn, ";",
                "samtools index", final.fn)
    if (verbose) {
        message("Preparing final output")
        message(cmd)
    }
    sys.res = system(cmd)
    if (sys.res) {
        file.remove(final.fn)
        stop("Error!")
    }

    return(final.fn)
}

#' @name loose.reads2
#' @title loose.reads2
#'
#' @description
#' 
#' Performs the following:
#' - loads reads from a specified window around loose end
#' - identifies any split reads and loads the "sides" of those reads as well
#' - loads the mates of reads around the loose end (by qname)
#' - realigns all of these reads
#'
#' The supplied BAM files are expected to have tag SA (indicated split alignments)
#' The functionality is identical to loose.reads; however, the BAMs are expected to be prefiltered to contain only the relevant QNAMEs via samtools to avoid slow filtering in R for reads with many mate windows
#'
#' @param tbam (character) path to tumor BAM file
#' @param taln (character) path realigned tumor bam
#' @param nbam (character) optional, path to normal BAM file, default=NULL
#' @param naln (character) optional, realigned normal bam
#' @param id (character)
#' @param filter (logical) return loose read pairs only? default FALSE
#' @param verbose optional, default=FALSE
#'
#' @return data.table of reads relaigned to the specified reference
#' @export
loose.reads2 = function(tbam, taln, nbam=NA, naln=NA, id="", filter=FALSE, verbose=FALSE){
    treads = .sample.spec2(tbam, verbose = verbose)
    realn = parse_realignment(treads, aln_bam = taln, filter = filter, verbose = verbose)
    realn$sample = id
    if(check_file(nbam) && check_file(naln)) {
        nreads = .sample.spec2(nbam, verbose = verbose)
        nrealn = parse_realignment(nreads, aln_bam = naln, filter = filter, verbose = verbose)
        nrealn$sample = paste0(id, "N")
        return(rbind(realn, nrealn, fill=TRUE, use.names=TRUE))
    }
    return(realn)
}

#' parse_realignment
#'
#' @description
#' takes two BAMs, one with all reads and mates, and one with realignments
#' and merges them, marking reads which did not realign uniquely to the reference
#'
#' @param reads (GRanges) reads from .sample_spec
#' @param aln_bam (character) bam realigned to reference
#' @param filter (logical) default FALSE
#' @param verbose (logical) default FALSE
#'
#' @return data.table of realigned reads
parse_realignment = function(reads, aln_bam, filter = FALSE, verbose = FALSE) {
    if (verbose) {
        message("Reading realigned reads from: ", aln_bam)
    }
    ## define standard chromosomes
    seqs = c(1:22, "X", "Y")
    bamseqs = seqnames(seqinfo(BamFile(aln_bam)))
    if(any(grepl("chr", bamseqs))) {
        seqs = gr.chr(seqs)
    }
    
    ## specifically we want alignment score here
    aln.grl = read.bam(bam = aln_bam, all = TRUE, pairs.grl = TRUE, tag="AS", isPaired = NA)
    aln.reads = as.data.table(unlist(aln.grl))

    ## trim /1 and /2 on qnames
    if (aln.reads[, .N]) {
        aln.reads[, qname := gsub("/[1|2]$", "", qname)]
        aln.reads[, R1 := bamflag(flag)[, "isFirstMateRead"]==1]
        aln.reads[, R2 := bamflag(flag)[, "isSecondMateRead"]==1]
    }
    
    ## recast reads as data table
    reads.dt = as.data.table(reads)
    ## setkeyv(reads.dt, c("qname", "R1"))
    ## identify which reads are unaligned
    if (verbose) {
        message("Identifying unaligned reads")
    }

    ## R1/R2 and qname should come from the original BAM,
    ## but the other characteristics of realn should come from single end BWA MEM

    ## reverse complement the negative strand seqs
    minus.strand.ix = which(bamflag(aln.reads[, flag])[, "isMinusStrand"]==1)
    aln.reads[, reading.frame := seq]
    aln.reads[minus.strand.ix, reading.frame := as.character(reverseComplement(DNAStringSet(reading.frame)))]

    ## for each aligned read, get the index of the original read and store as query.id
    qid = match(paste(aln.reads$qname, aln.reads$reading.frame), paste(reads.dt$qname, reads.dt$seq))
    ## paste(reads.dt$qname, reads.dt$R1, sep = "_"))
    aln.reads = aln.reads[, query.id := qid][!is.na(query.id)]
    aln.reads[(!is.na(query.id)), ":="(R1 = reads.dt$R1[query.id],
                                       R2 = reads.dt$R2[query.id],
                                       flag = reads.dt$flag[query.id])] ## keep original flag?

    
    ## for unaligned reads, also get the index and store that as query.id
    unaln.reads.ix = which(is.na(match(reads.dt$seq, aln.reads$reading.frame)))
    ##reads.dt[aln.reads[, .(qname, R1)], nomatch = NA][, which(is.na(seq))]

    unaln.reads = reads.dt[unaln.reads.ix][, ":="(seqnames = "*",
                                               start = 1,
                                               end = 0,
                                               flag = ifelse(R1, 69, 113),
                                               mapq = 0,
                                               query.id = unaln.reads.ix)]
    
    ## return concatenated aligned reads and unaligned reads
    realn = rbind(aln.reads, unaln.reads, fill = TRUE, use.names = TRUE)

    ## MAPQ is zero if not on a standard chromosome
    realn$mapq = as.integer(realn$mapq)
    realn[, mapq := ifelse(!(seqnames %in% seqs), as.integer(0), mapq), by=query.id]
    realn = realn[rev(order(mapq))][!duplicated(query.id), ]
    realn[, MQ := ifelse(rep(.N, .N)==1, as.integer(NA), c(mapq[-1], mapq[1])), by=qname]
    cols = c(colnames(realn)[colnames(realn) %in% colnames(as.data.table(reads))], "reading.frame", "AS")
    realn = realn[, cols, with=F]

    ## annotate whether the read belongs to a loose read pair
    lqn = realn[mapq>50 & (is.na(MQ) | MQ < 1), qname]
    if(filter){
        realn = realn[qname %in% lqn,]
        realn[, loose.pair := TRUE]
    } else {
        realn[, loose.pair := qname %in% lqn]
    }
    ## mapq is the probability of given read aligning
    ## MQ is mapq of the mate
    realn[, high.mate := mapq>50 & (is.na(MQ) | MQ < 1)]
    return(realn)
}
