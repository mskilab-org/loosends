#' @importFrom bamUtils read.bam

#' @name bowtie_aln
#' @title bowtie_aln
#'
#' @description
#'
#' Perform unpaired (single end) realignment of reads to either a supplied sequence or reference FASTA
#' 
#' @param reads (data.table with columns qname, seq)
#' @param ref.seq (data.table with columns qname, seq)
#' @param ref.dir
#' @param ref.basename
#' @param outdir
#' @param verbose
#' 
#' @return GRanges with alignment
bowtie_aln = function(reads,
                      ref.seq = NULL,
                      ref.dir = NULL,
                      ref.basename = NULL,
                      outdir = "./",
                      verbose = FALSE) {

    ## if there are no reads just return empty GRanges
    if (!reads[, .N]) {
        message("No reads!")
        return(GRanges())
    }

    ## check if we need to create/index FASTA
    if (!is.null(ref.seq)) {
        ref.dir = outdir
        ref.basename = "target"
        fasta.fn = paste0(outdir, "/", ref.basename, ".fasta")
        if (verbose) { message("Writing reference to .fasta: ", fasta.fn) }
        fasta.fn = seqs_to_fasta(ref.seq, fasta.fn, verbose = verbose)
        if (verbose) { message("Indexing fasta") }
        ref.dir = bowtie_index(fasta.fn = fasta.fn, ref.id = ref.basename, verbose = verbose)
    } else {
        if (is.null(ref.dir) | is.null(ref.basename)) {
            stop("ref.dir and ref.basename must be supplied!")
        }
        if (verbose) {
            message("Using reference: ", ref.basename)
            message("Located in directory: ", ref.dir)
        }
    }

    if (verbose) { message("Writing reads to fastq") }
    reads.fastq.fn = seqs_to_fasta(reads[, .(qname, seq = reading.frame)],
                                   fasta.fn = paste0(outdir, "/reads.fastq"),
                                   format = "fastq",
                                   verbose = verbose)
    
    ## perform alignment
    if (verbose) {message("Starting alignment")}
    
    currdir = normalizePath(getwd())
    if (verbose) {
        message("Current directory: ", currdir)
        message("Switching to new directory: ", ref.dir)
    }
    setwd(ref.dir)

    aligned.fn = paste0(normalizePath(outdir), "/aln.sam")
    
    cmd = paste("module load bowtie2;",
                "bowtie2 --very-sensitive-local -x", ref.basename,
                "-U", reads.fastq.fn,
                "-S", aligned.fn)

    if (verbose) {
        message("RUNNING: ", cmd)
    }

    sys.res = system(cmd)
    if (sys.res) {
        file.remove(aligned.fn)
        stop("Error!")
    }

    if (verbose) {
        message("Going back to OG directory: ", currdir)
    }
    setwd(currdir)

    ## sort
    final.fn = normalizePath(paste0(currdir, "/aln.bam"))
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

    ## read BAM file and replace qnames
    aln.grl = read.bam(bam = final.fn, all = TRUE, pairs.grl = TRUE, tag="AS", isPaired = NA)
    aln.reads = as.data.table(unlist(aln.grl))

    if (!nrow(aln.reads)) {
        return(GRanges())
    }

    ## remove unmapped reads
    aln.reads = aln.reads[!(start == 1 & end == 0),]

    if (!nrow(aln.reads)) {
        return(GRanges())
    }

    aln.reads[, ix := as.numeric(as.character(qname))]
    aln.reads[, qname := reads[ix, qname]]
    aln.reads[, R1 := reads[ix, R1]]

    ## check if seqnames should also be replaced
    if (!is.null(ref.seq)) {
        aln.reads[, seqnames.ix := as.numeric(as.character(seqnames))]
        aln.reads[, seqnames := as.character(ref.seq[seqnames.ix, qname])]
    }

    return(dt2gr(aln.reads))
}

#' @name bowtie_index
#' @title bowtie_index
#'
#' @param fasta.fn (character) path to fasta
#' @param ref.id (character) basename of reference, default "target"
#' @param verbose (logical) default FALSE
#'
#' @return directory containing indexed FASTA
bowtie_index = function(fasta.fn = NA_character_, ref.id = "target", verbose = FALSE) {

    if (!check_file(fasta.fn)) {
        stop("fasta does not exist: ", fasta.fn)
    }

    currdir = normalizePath(getwd())
    newdir = normalizePath(dirname(fasta.fn))
    if (verbose) {
        message("Current directory: ", currdir)
        message("Switching to new directory for indexing: ", newdir)
    }
    setwd(newdir)

    ## cmd = paste("bowtie2-build", basename(fasta.fn), ref.id)
    cmd = paste("module load bowtie2; bowtie2-build", basename(fasta.fn), ref.id)
    if (verbose) {
        message("RUNNING: ", cmd)
    }

    sys.res = system(cmd)

    if (sys.res) {
        stop("Error!")
    }

    if (verbose) {
        message("Finished indexing reference!")
        message("Switching back to directory: ", currdir)
    }

    return(newdir)
}

#' @name seqs_to_fasta
#' @title seqs_to_fasta
#'
#' @param seqs.dt (data.table) must have columns qname and fastq
#' @param fasta.fn (character) fasta file name to write to
#' @param format (fasta or fastq)
#' @param verbose (logical) default FALSE
#'
#' @return fasta.fn
seqs_to_fasta = function(seqs.dt, fasta.fn, format = "fasta", verbose = FALSE) {

    outdir = normalizePath(dirname(fasta.fn))
    if (!dir.exists(outdir)) {
        if (verbose) { message("Creating output directory: ", outdir) }
        dir.create(outdir, recursive = TRUE)
    }

    if (!("seq" %in% colnames(seqs.dt))) {
        stop("seqs.dt must contain column $seq")
    }

    if (!("qname" %in% colnames(seqs.dt))) {
        stop("seqs.dt must contain column $qname")
    }

    ## use the number as qname
    if (verbose) {
        message("Creating DNAStringset and writing to file: ", normalizePath(fasta.fn))
    }
    seq.vec = setNames(object = seqs.dt[, seq], nm = as.character(1:seqs.dt[, .N]))
    ##seq.vec = setNames(object = seqs.dt[, seq], nm = seqs.dt[, qname])##as.character(1:seqs.dt[, .N]))
    seq.strings = Biostrings::DNAStringSet(seq.vec)

    Biostrings::writeXStringSet(x = seq.strings, filepath = normalizePath(fasta.fn), format = format)

    return(normalizePath(fasta.fn))
}
