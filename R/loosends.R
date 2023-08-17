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
#' @import readsupport
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import stringr

#' @name prep_loose_ends
#' @title prep_loose_ends
#'
#' @param li (GRanges, character, or data.table)
#' @param id (character) sample ID
#' 
#' @return data.table of loose ends with columns sample and leix
prep_loose_ends = function(li, id = "") {
    if (inherits(li, "character")) {
        le.dt = as.data.table(parse.gr(li))
    } else if (inherits(li, "GRanges")) {
        le.dt = as.data.table(li)
    } else if (inherits(li, "data.table")) {
        le.dt = copy(li)
    }

    if (is.null(le.dt$sample)) {
        le.dt[, sample := id]
    }
    if (is.null(le.dt$leix)) {
        le.dt[, leix := 1:.N]
    }

    return(le.dt)
}

    
#' @name prep_loose_reads
#' @title prep_loose_reads
#'
#' Get loose reads data table ready for assembly
#'
#' @param li (data.table from prep_loose_ends)
#' @param loose.reads.dt (from loosereads_wrapper)
#'
#' @return data.table
#' @export
prep_loose_reads = function(li, loose.reads.dt) {
    if (nrow(li) > 1) {
        stop("More than 1 loose end provided; did you mean 'process.loose.ends'?")
    }
    if (nrow(li)==0) {
        stop("0 loose ends provided")
    }

    if (is.null(li$sample)) {
        stop("li must contain column $sample")
    }

    if (is.null(li$leix)) {
        stop("li must contain column $leix")
    }

    if (!nrow(loose.reads.dt)) {
        return(loose.reads.dt)
    }

    if (is.null(loose.reads.dt$sample)) {
        stop("loose.reads.dt must contain column $sample")
    }

    ## get loose ends objects ready
    ri = copy(loose.reads.dt)
    ri$sample = as.character(ri$sample)

    ## add leix to both ri and li because this is needed for saving samples
    ri$leix = li$leix

    ## denote sample vs. control
    ri[, track := paste(ifelse(sample == li$sample, "sample", "control"),
                        ifelse(strand == "+", "for", "rev"), sep=".")]
    ri[, concord := !(loose.pair) & .N == 2 & length(unique(seqnames)) == 1 & strand[R1] != strand[R2] & strand[start == min(start)]=="+" & min(start) + 3e3 > max(start), by=qname]
    ri[, anchor := (loose.pair & high.mate) | ( !(loose.pair) & mapq > 50 & !(concord))]
    ri$seq = as.character(ri$seq)
    ri$start = as.integer(ri$start)
    ri$end = as.integer(ri$end)
    ri$flag = as.integer(ri$flag)

    ## remove short sequences
    ## median.nchar = median(nchar(ri$seq))
    ## ri = ri[(nchar(seq) >= median.nchar),]

    ## do this regardless
    ri$reading.frame = ri$seq
    ri[strand == "-", reading.frame := as.character(reverseComplement(DNAStringSet(reading.frame)))]

    return(ri)
}

#' @name minimap_index
#' @title minimap_index
#'
#' @param fasta (character) path to fasta to pre-index
#' @param outdir (character) path to target directory (created if not already exists)
#' @param verbose (logical) default FALSE
#'
#' @return path to minimap indexed reference (.mmi)
minimap_index = function(fasta = NA_character_, outdir = "./", verbose = FALSE) {
    if (!dir.exists(outdir)) {
        if (verbose) message("Creating output directory")
        dir.create(outdir, recursive = TRUE)
    }

    if (!check_file(fasta)) {
        stop("Fasta does not exist: ", fasta)
    }

    target.fn = paste0(outdir, "/", "target.mmi")
    cmd = paste("minimap2 -x sr -d", target.fn, fasta)
    sys.res = system(cmd)

    if (sys.res) {
        file.remove(target.fn)
        stop("Error!")
    }

    return(target.fn)
}
