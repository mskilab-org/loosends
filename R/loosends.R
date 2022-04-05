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

## #' @name call_loose_end_wrapper
## #' @title call_loose_end_wrapper
## #'
## #' @description
## #'
## #' call loose end given reads
## #'
## #' @param le.dt (from prep_loose_ends)
## #' @param reads.dt (from loosereads)
## #' @param concat.bwa (BWA from RSeqLib) concatenated fasta
## #' @param human.bwa (BWA from RSeqLib)
## #' @param pad (window for local assembly)
## #' @param mix.tn (logical) mix tumor and normal reads? default TRUE
## #' @param minimap (logical) default FALSE
## #' @param max.reads (numeric) max reads per loose end before downsampling (default 5000)
## #' @param verbose
## #'
## #' @return list with three items:
## #' $call table with calls
## #' $filtered.contigs data.table with filtered contigs
## #' $all.contigs data.table with all contigs
## #'
## #' @export
## call_loose_end_wrapper = function(id = "",
##                                   le.dt = data.table(),
##                                   reads.dt = data.table(),
##                                   concat.bwa = NULL,
##                                   human.bwa = NULL,
##                                   concat.fn = NA_character_,
##                                   pad = 5e3,
##                                   mix.tn = TRUE,
##                                   minimap = FALSE,
##                                   outdir = "./",
##                                   max.reads = 5000,
##                                   verbose = FALSE) {

##     if (is.null(concat.bwa)) { stop("Must supply BWA object as concat.bwa") }
##     if (is.null(human.bwa)) { stop("Must supply BWA object as human.bwa") }
        
##     empty.res = list(call = data.table(), filtered.contigs = data.table(), all.contigs = data.table())

##     le.dt = prep_loose_ends(li = le.dt, id = id)

##     if (!le.dt[, .N]) {
##         return(empty.res)
##     }

##     if (!reads.dt[, .N]) {
##         return(empty.res)
##     }

##     full.res = lapply(1:le.dt[, .N],
##                       function(ix) {
##                           sel = dt2gr(reads.dt) %^% (dt2gr(le.dt[ix,]) + pad)
##                           if (any(sel)) {
##                               qns = reads.dt[sel, qname]
##                               this.reads.dt = reads.dt[qname %in% qns,]
##                           } else {
##                               this.reads.dt = reads.dt
##                           }
##                           ri = prep_loose_reads(li = le.dt[ix,],
##                                                 loose.reads.dt = this.reads.dt)

##                           sub.res = call_loose_end(li = le.dt[ix,],
##                                                    ri = ri,
##                                                    concat.bwa = concat.bwa,
##                                                    human.bwa = human.bwa,
##                                                    pad = pad,
##                                                    mix.tn = mix.tn,
##                                                    concat.fn = concat.fn,
##                                                    minimap = minimap,
##                                                    outdir = paste0(outdir, "/", ix),
##                                                    max.reads = max.reads,
##                                                    verbose = verbose)
##                           return(sub.res)
##                       })

##     full.call = lapply(full.res, function(x) {x$call}) %>%
##         rbindlist(fill = TRUE, use.names = TRUE)
##     full.filtered = lapply(full.res, function(x) {x$filtered.contigs}) %>%
##         rbindlist(fill = TRUE, use.names = TRUE)
##     full.contigs = lapply(full.res, function(x) {x$all.contigs}) %>%
##         rbindlist(fill = TRUE, use.names = TRUE)

##     return(list(call = full.call, filtered.contigs = full.filtered, all.contigs = full.contigs))
## }

## #' @name build.from.win
## #'
## #' assembles strand-specific reads and their mates into contigs and counts read support per contig within a single seed window
## #' @param win GRanges seed window
## #' @param ri data.table reads and their mates
## #' @param tracks optional, character which tracks (sample + strand) to assemble, default all
## #' @param align.thres (numeric) width of alignment / number of non-N characters
## build.from.win = function(win, ri, tracks=NULL, verbose = FALSE, align.thres = 0.9){
##     if(!("track" %in% colnames(ri))){
##         ri[, track := ifelse(strand=="+", "for", "rev")]
##     }
##     if(is.null(tracks)) tracks = ri[, unique(track)]

##     if (verbose) {
##         message("Tracks used for assembly: ")
##         message(paste(tracks, collapse = '\n'))
##     }
##     ri[track %in% tracks, {
##         t = track[1]
##         ri[, seed := dt2gr(ri) %N% win > 0 & track==t]
##         rtmp = ri[qname %in% qname[seed]]
##         rtmp[(seed)][is.dup(qname), seed := 1:.N == 1, by=qname] ## in case of a foldback where both reads are "seed", choose one as the anchor
##         rtmp[, seed.frame := ifelse(seed, reading.frame, as.character(reverseComplement(DNAStringSet(reading.frame))))]
##         if(nrow(rtmp) < 5){
##             data.table(peak = gr.string(win[,c()]), seq = as.character(NA), good.assembly=FALSE, cov=ri[(seed), .N], mapq60=ri[(seed) & mapq==60, .N], unassembled.reads=nrow(rtmp))
##         } else{
##             srf = rtmp[, seed.frame]
            
##             if (verbose) {
##                 message("Starting assembly with ", length(srf), " reads")
##             }
##             f = Fermi(srf, assemble=T)
##             contigsf = contigs(f)
##             unassembled.reads = length(srf)

##             if (verbose) {
##                 message("Initial number of contigs: ", length(contigsf))
##             }
            
##             if (length(contigsf)) {
##                 reassemble = TRUE
##                 while(reassemble){

##                     ## using current contigs, count the number of unassembled reads
##                     ctigs = BWA(seq=contigsf)
##                     aln = ctigs[srf]

##                     ## qnames not included in alignment
##                     s1 = which(!(seq_along(srf) %in% as.integer(aln$qname)))

##                     ## qnames with "bad" alignments (e.g. alignment width < 0.9 * #non-N characters)
##                     aln.dt = as.data.table(aln)
##                     align.frac.dt = aln.dt[, .(best.aln = max(width / nchar(gsub("N", "", seq)))), by = qname]
##                     s2 = align.frac.dt[best.aln < align.thres, as.integer(qname)]
##                     ## s2 = as.data.table(aln)[, max(width / nchar(gsub("N", "", seq))), by=qname][V1 < 0.9, as.integer(qname)]
##                     ## if the number of unassembled reads does not change from the previous iteration...
##                     if(length(s1) + length(s2) == unassembled.reads | length(s1) + length(s2) < 7) {
##                         reassemble = FALSE
##                     } else {
##                         unassembled.reads = length(s1) + length(s2)

##                         if (verbose) {
##                             message("Reassembling. Number of unassembled reads: ", unassembled.reads)
##                         }
##                         ## if(unassembled.reads > 7){
##                         ## message("Starting assembly with ", unassembled.reads, " reads")
##                         c2 = Fermi(srf[c(s1, s2)], assemble=T)
##                         if(length(contigs(c2))){
##                             contigsf = c(contigsf, contigs(c2))
##                         } else {
##                             reassemble = FALSE
##                         }
##                     }
##                 }
##                 goodassembly = rep(TRUE, length(contigsf))
##                 for(i in seq_along(contigsf)){
##                     ctig = BWA(seq = contigsf[i])
##                     ra = ctig[srf]
##                     sc = table(factor(strand(ra), c("+", "-")))
##                     if(!any(sc[1:2] == 0)) goodassembly[i] = FALSE
##                     if(sc["+"] < sc["-"]) contigsf[i] = as.character(reverseComplement(DNAStringSet(contigsf[i])))
##                 }
##                 data.table(peak = gr.string(win[,c()]), seq = contigsf, good.assembly=goodassembly, cov=ri[(seed), .N], mapq60=ri[(seed) & mapq==60, .N], unassembled.reads=unassembled.reads)
##             } else{
##                 data.table(peak = gr.string(win[,c()]), seq = as.character(NA), good.assembly=FALSE, cov=ri[(seed), .N], mapq60=ri[(seed) & mapq==60, .N], unassembled.reads = length(srf))
##             }
##         }
##     }, by=track]
## }

## #' match.seq
## #'
## #' returns logical vector of length subject indicating whether matches were found with any query
## #' @param query PDict of query sequences to match
## #' @param subject DNAStringSet of sequences to parse for matches

## match.seq = function(query, subject)
## {
##     if (is.null(names(subject)))
##         names(subject) = 1:length(subject)

##     if (any(duplicated(names(subject))))
##     {
##         warning('Names of subject sequences have duplicates, deduping')
##         names(subject) = dedup(names(subject))
##     }

##     totals = base::lengths(Biostrings::vwhichPDict(query, subject)) > 0
##     return(totals)
## }

#' @name eighteenmer
#' @title eighteenmer
#'
#' @description
#' 
#' generate a PDict of telomeric 18mers
#' 
#' @param gorc 'g' returns 18mers of G-heavy strand only 'c' returns 18mers of C-heavy strand only 'both' returns 18mers of G-heavy strand and 18mers of C-heavy strand, default='both'
eighteenmer = function(gorc=c('g', 'c', 'both')){
    if(any(!(gorc %in% c('g', 'c', 'both')))) stop("Allowed values of gorc are 'g','c','both'")
    if(length(gorc)>1) gorc='both'
    pattern = system.file('extdata', 'telomeres.fa', package='loosends')
    query = tryCatch(readDNAStringSet(pattern), error = function(e) NULL)
    if (is.null(query))
    {
        query = DNAStringSet(readLines(pattern))
    }
    motifs = unname(as.character(query))
    if(gorc!='c'){
        g.motifs = motifs[grepl("GGG", motifs)]
        gquery = as.data.table(expand.grid(i=g.motifs,j=g.motifs,k=g.motifs,l=g.motifs,s=1:6))[, seq := substr(paste0(i, j, k,l), s, s+17)][, unique(seq)]
        gqueries = PDict(gquery)
        if(gorc=='g') return(gqueries)
    }
    if(gorc!='g'){
        c.motifs = motifs[grepl("CCC", motifs)]
        cquery = as.data.table(expand.grid(i=c.motifs,j=c.motifs,k=c.motifs,l=c.motifs,s=1:6))[, seq := substr(paste0(i, j, k,l), s, s+17)][, unique(seq)]
        cqueries = PDict(cquery)
        if(gorc=='c') return(cqueries)
    }
    return(PDict(c(gquery, cquery)))
}

#' @name munch
#' @title munch
#'
#' @description
#' 
#' returns logical vector of length reads
#' indicating does read contain any match to query
#' searches for matches on same strand as input seq
#' 
#' @param reads GRanges or data.table containing $seq field with character sequence
#' @param query optional, PDict of motifs, default= 18mers of telomeric motifs (both strands)
munch = function(reads, query=NULL) {
    if(is.null(query)) query = eighteenmer()
    seq = reads$seq
    seq_frame_ref = DNAStringSet(seq)
    return(match.seq(query, seq_frame_ref))
}

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
    

## #' @name grab_ref_obj
## #' @title grab_ref_obj
## #'
## #' @param ref.dir (character) path to reference directory
## #'
## #' @return list with names rep, human, polyA, microbe
## grab_ref_obj = function(ref.dir = NA_character_) {

##     if (!dir.exists(ref.dir)) {
##         stop("Supplied directory does not exist")
##     }

##     human = BWA(fasta=paste(ref.dir, "human_g1k_v37_decoy.fasta", sep="/"))
##     rep = BWA(fasta=paste(ref.dir, "mskilab_combined_TraFicv8-3_satellites.fa", sep="/"))
##     polyA = BWA(fasta=paste(ref.dir, "PolyA.fa", sep="/"))
##     microbe = BWA(fasta=paste(ref.dir, "human_g1k_v37.withviral.fasta", sep="/"),
##                   keep_sec_with_frac_of_primary_score=0.2)

##     ref.obj = list(human=human, rep=rep, polyA=polyA, microbe=microbe)
##     return(ref.obj)
## }

## #' @name call_loose_end
## #' @title call_loose_end
## #'
## #' @description
## #'
## #' This function calls a loose end given reads from near to the loose end
## #'
## #' @param li (data.table) data.table coercible to GRanges, must have metatdata $leix and $sample
## #' @param ri (data.table) loose reads table prepped (e.g. from prep_loose_reads)
## #' @param concat.bwa (BWA object from RSeqLib)
## #' @param human.bwa (BWA object from RSeqLib)
## #' @param pad (numeric) window size for local assembly, default 1 kbp
## #' @param mix.tn (logical) mix tumor/normal reads before assembly
## #' @param max.reads (logical) maximum number of reads before downsampling (default 5000)
## #' @param verbose
## #'
## #' @return list with entries
## #' - all.contigs (data.table)
## #' - wide.contigs (data.table or NULL)
## #' - call (data.table)
## call_loose_end = function(li, ri,
##                           concat.bwa = NULL, ## BWA object with concatenated reference
##                           human.bwa = NULL, ## BWA object with just human reference
##                           concat.fn = "~/git/loosends/inst/extdata/hg19_loosends/concatenated_references_deduped.fasta",
##                           pad = 5e3,
##                           mix.tn = TRUE,
##                           max.reads = 2500,
##                           outdir = "./",
##                           minimap = FALSE,
##                           verbose = FALSE) {
    
##     ## human = ref_obj$human
##     ## rep = ref_obj$rep
##     ## polyA = ref_obj$polyA
##     ## microbe = ref_obj$microbe
    
##     uannot = readRDS(system.file('extdata', '101.unmappable.annotations.rds', package='loosends'))

##     id = li$sample

##     if (!dir.exists(outdir)) {
##         if (verbose) message("Creating output directory")
##         dir.create(outdir, recursive = TRUE)
##     }

##     if (verbose) {
##         message("Generating call")
##     }

##     .build.tigs = function(ri, pp, id, leix, verbose = TRUE) {
##         if(verbose) message("assembling contigs")
##         out.dt = rbindlist(lapply(1:length(pp), function(i){
##             win = pp[i]
##             build.from.win(win, ri, verbose = verbose)
##         }), use.names=TRUE, fill=TRUE)
##         if("seed" %in% colnames(ri)) ri$seed = NULL
##         out.dt$sample = rep(id, nrow(out.dt))
##         out.dt$leix = rep(leix, nrow(out.dt))
##         out.dt$somatic = rep(somatic, nrow(out.dt))
##         out.dt = out.dt[!is.na(seq)]

##         if(verbose) message(paste("assembled", nrow(out.dt), "contigs across all tracks"))
##         if(nrow(out.dt)){
##             out.dt[, qname := 1:.N, by=.(track, sample, leix)]
##             out.dt[, qname := paste(sample, leix, track, qname, sep="_")]
##         } else{
##             out.dt$qname = character()
##         }

##         ## save contigs as FASTA
##         if (out.dt[, .N]) {
##             unique.tigs = unique(out.dt, by = "qname")
##             tigs.vec = setNames(object = unique.tigs[, seq], nm = unique.tigs[, qname])
##             tigs.strings = Biostrings::DNAStringSet(tigs.vec)
##             Biostrings::writeXStringSet(tigs.strings, paste0(outdir, "/", "contigs.fasta"))
##         }

##         if (verbose) {
##             message("Aligning contigs to concatenated reference")
##         }

##         if (minimap) {
##             if (verbose) message("Using minimap for contig alignment")

##             contigs.aln.bam.fn = paste0(outdir, "/", "contigs.aln.bam")
##             contigs.aln.sorted.bam.fn = paste0(outdir, "/", "contigs.aln.sorted.bam")
##             contigs.aln.sam.fn = paste0(outdir, "/", "contigs.aln.sam")

##             if (verbose) message("starting short read single end alignment")
##             cmd = paste("minimap2 -ax sr --secondary=yes",
##                         concat.fn,
##                         paste0(outdir, "/", "contigs.fasta"),
##                         ">",
##                         contigs.aln.sam.fn)

##             sys.res = system(cmd)
##             if (sys.res) {
##                 file.remove(contigs.aln.sam.fn)
##                 stop("Error!")
##             }

##             cmd = paste("samtools view -Sb", contigs.aln.sam.fn, ">", contigs.aln.bam.fn)
##             sys.res = system(cmd)
            
##             if (sys.res) {
##                 file.remove(contigs.aln.bam.fn)
##                 stop("Error!")
##             }

##             cmd = paste("samtools sort", contigs.aln.bam.fn, ">", contigs.aln.sorted.bam.fn, ";",
##                         "samtools index", contigs.aln.sorted.bam.fn)
##             sys.res = system(cmd)
            
##             if (sys.res) {
##                 file.remove(contigs.aln.bam.sorted.fn)
##                 stop("Error!")
##             }

##             ## read BAM from single end alignment
##             aln.gr = read.bam(bam = contigs.aln.sorted.bam.fn,
##                               all = TRUE,
##                               pairs.grl = FALSE,
##                               tag="AS",
##                               isPaired = NA)
##             ralns.og = aln.gr
##             ## don't include unaligned reads?
##             if (length(ralns.og)) {
##                 ralns.og = ralns.og %Q% (!is.na(qname))
##             }

##             ## garbage collect afterwards to avoid memory issues??
##             ## hack hack
##             gc()
##         } else {
##             ralns.og = concat.bwa[out.dt$seq]
##             if (length(ralns.og)) {
##                 values(ralns.og)[, "qname"] = out.dt[as.integer(ralns.og$qname), qname]
##             }
##         }

##         ## TODO: store this as a package global variable
##         keepseq = seqnames(seqinfo(concat.bwa))[c(1:87, 148:6398)]
##         rep.seqs = seqnames(seqinfo(concat.bwa))[1:60]
##         poly.seqs = seqnames(seqinfo(concat.bwa))[61]
##         human.seqs = seqnames(seqinfo(concat.bwa))[62:86]
##         viral.seqs = seqnames(seqinfo(concat.bwa))[148:6398]

        
##         keep.cols = c("cigar", "flag", "mapq", "AS", "qname")

##         if (length(ralns.og)) {
##             ralns.dt = as.data.table(ralns.og[, keep.cols]  %Q% (as.character(seqnames(ralns.og)) %in% keepseq))

##             ## transfer original values based on qname
##             setnames(ralns.dt, old = "qname", new = "qname.og")
##             out.dt.cols = setdiff(colnames(out.dt), colnames(ralns.dt))
##             ralns.dt = cbind(ralns.dt, out.dt[match(ralns.dt$qname.og, qname), ..out.dt.cols])

##             ralns.dt[as.character(seqnames) %in% rep.seqs, c_type := "rep"]
##             ralns.dt[as.character(seqnames) %in% poly.seqs, c_type := "polyA"]
##             ralns.dt[as.character(seqnames) %in% human.seqs, c_type := "human"]
##             ralns.dt[as.character(seqnames) %in% viral.seqs, c_type := "viral"]
##             ralns.dt[, c_spec := c_type]
##             ralns.dt[c_type == "rep", c_spec := dunlist(strsplit(as.character(seqnames), "#", fixed=T))[rev(!duplicated(rev(listid)))][order(as.integer(listid)), V1]]
##             calns = ralns.dt[!(is.na(cigar) | is.na(AS) | is.na(flag) | is.na(mapq))]
##         } else {
##             calns = as.data.table(ralns.og)
##         }

##         ## if(verbose) message("aligning contigs to 4 references")
##         ## ureps = rep[out.dt$seq]
##         ## preps = polyA[out.dt$seq]
##         ## mreps = microbe[out.dt$seq]
##         ## hreps = human[out.dt$seq]
##         ## virals = seqlevels(microbe)[85:length(seqlevels(microbe))]
        
##         ## if(!(length(ureps))){
##         ##     ureps$cigar = character()
##         ##     ureps$mapq = character()
##         ##     ureps$AS = integer()
##         ##     ureps$flag = character()
##         ## }
##         ## if(!(length(preps))){
##         ##     preps$cigar = character()
##         ##     preps$mapq = character()
##         ##     preps$AS = integer()
##         ##     preps$flag = character()
##         ## }
##         ## if(!(length(hreps))){
##         ##     hreps$cigar = character()
##         ##     hreps$mapq = character()
##         ##     hreps$AS = integer()
##         ##     hreps$flag = character()
##         ## }
##         ## if(!(length(mreps))){
##         ##     mreps$cigar = character()
##         ##     mreps$mapq = character()
##         ##     mreps$AS = integer()
##         ##     mreps$flag = character()
##         ## }
##         ## keep.cols = c("cigar", "flag", "mapq", "AS")
##         ## values = rbind(out.dt[as.integer(ureps$qname)], out.dt[as.integer(preps$qname)], out.dt[as.integer(hreps$qname)], fill=T, use.names=T)
##         ## ralns = rbind(as.data.table(ureps[, keep.cols]), as.data.table(preps[, keep.cols]), as.data.table(hreps[, keep.cols]))
##         ## ralns = cbind(ralns, values)
##         ## if(nrow(out.dt)){
##         ##     rch = cgChain(ralns)
##         ##     good.ids = c(as.character(seqnames(gaps(rch$x) %Q% (strand=="+"))), setdiff(out.dt$qname, ralns$qname))
##         ##     mreps$query.id  = as.integer(mreps$qname)
##         ##     mreps$qname = out.dt[mreps$query.id, qname]
##         ##     vreps = mreps %Q% (seqnames %in% virals) %Q% (qname %in% good.ids)
##         ##     if(verbose & length(vreps)) message("adding viral alignments")
##         ##     valns = cbind(as.data.table(vreps[, keep.cols]), out.dt[vreps$query.id])
##         ##     calns = rbind(ralns, valns)
##         ##     calns[, c_type := c(rep("rep", length(ureps)), rep("polyA", length(preps)), rep("human", length(hreps)), rep("viral", length(vreps)))]
##         ##     calns[, c_spec := c_type]
##         ##     calns[c_type == "rep", c_spec := dunlist(strsplit(as.character(seqnames), "#", fixed=T))[rev(!duplicated(rev(listid)))][order(as.integer(listid)), V1]]
##         ## } else{
##         ##     calns = ralns
##         ## }
##         calns$somatic = rep("somatic", nrow(calns))
##         calns$mapq = as.integer(calns$mapq)

##         ## fill in dummy values
##         calns[, ":="(read.cov = NA, map60.cov = NA, tumor.frac = NA,
##                      tumor.support = NA, normal.support = NA, total.support = NA)]


##         if(verbose) message("parsing contigs for telomeric matches")
##         all.contigs = copy(calns)
##         if (nrow(all.contigs)) {
##             all.contigs[c_type == "human", c_spec := ifelse(seqnames %in% c(1:22, "X", "Y"), "human", "unassembled")]
##             all.contigs$cmer = munch(all.contigs, eighteenmer('c'))
##             all.contigs$gmer = munch(all.contigs, eighteenmer('g'))

##             if(verbose) message("parsing alignments for unmappable repeat overlaps")
##             coord.calls = copy(all.contigs[c_type=="human"])[, !c("c_type", "c_spec")]
##             ov = gr.findoverlaps(dt2gr(coord.calls), uannot)
##             if(length(ov)){
##                 subj = ov$subject.id
##                 strand(ov) = coord.calls[ov$query.id, strand]
##                 values(ov) = values(dt2gr(coord.calls)[ov$query.id])
##                 ov$c_type = "rep"
##                 ov$c_spec = uannot[subj]$c_spec
##                 coord.calls = as.data.table(ov)
##                 return(rbind(all.contigs, coord.calls, fill=T, use.names=T))
##             }
##         }
##         return(all.contigs)
##     }

##     if (verbose) {
##         message("Building contigs")
##     }

##     somatic = as.logical(nrow(ri[grepl("control", track)]))
##     wholseed = dt2gr(li)[,c()] + pad ##dt2gr(li)[,c()]+1e3
##     if (verbose) {
##         message("Using window size: ", pad, " bp")
##     }
##     pp = gr.stripstrand(gr.tile(wholseed, 200) %Q% (width == 200))

##     ## check if contigs should be assembled mixing tumor and normal reads
##     ri.input = copy(ri)
##     if (mix.tn) {
##         if (verbose) {
##             message("Mixing tumor/normal reads before assembly")
##         }
##         ri.input$track = NULL ## set track to null so that assembly is sample-agnostic
##     }
##     ## check if reads need to be downsampled before assembly
##     if (ri.input[, .N]) {
##         n.qn = length(unique(ri.input[, qname]))
##         if (verbose) {
##             message("Number of unique read qnames for assembly: ", n.qn)
##         }
##         if (n.qn > max.reads) {
##             if (verbose) {
##                 message("Downsampling reads to max # qnames: ", max.reads)
##             }
##             keep.qn = sample(unique(ri.input[, qname]), max.reads, replace = FALSE)
##             ri.input = ri.input[qname %in% keep.qn,]
##         }
##     }
##     all.contigs = .build.tigs(ri.input, pp, li$sample, li$leix, verbose = verbose)

##     ## also save the filtered contigs
##     if (verbose) {
##         message("Generating call")
##     }

##     if (mix.tn) {
##         res = caller(li, all.contigs, human = human.bwa, 
##                      return.contigs = TRUE, pad = 5e3)
##     } else {
##         ## if not mixed, call should only be generated from tumor contigs
##         res = caller(li, all.contigs[grepl('sample', track),],
##                      human = human.bwa,
##                      return.contigs = TRUE)
##     }


##     ## add junction annotation to cal
##     call = res$call
##     filtered.contigs = res$contigs
##     ## HACK HACK
##     bowtie = max(ri[, mapq], na.rm = TRUE) < 60 ## you know it's bowtie if mapqs only go up to 44
##     if (bowtie) {
##         recall.res = read.based(li, ri, pad = pad, human = human.bwa, return.contigs = TRUE, mapq.thresh = 30)
##     } else {
##         recall.res = read.based(li, ri, pad = pad, human = human.bwa, return.contigs = TRUE)
##     }
##     disc = recall.res$call
##     recall = (!call$missedj & disc$missedj) | (!call$complex & disc$complex)
##     call[, missedj := missedj | disc$missedj]
##     call[, complex := complex | disc$complex]
##     call[junction == "" & disc$junction != "", junction := disc$junction]
##     call[junction != "" & disc$junction != "", junction := {
##         cgr = parse.gr(junction)
##         dgr = parse.gr(disc$junction)
##         s = gr.reduce(cgr[1], dgr[1], ignore.strand=F)
##         m = gr.reduce(cgr[-1], dgr[-1], ignore.strand=F)
##         m = gr.start(GenomicRanges::reduce(m + 200) - 200, ignore.strand=F)
##         paste(c(gr.string(s), gr.string(m)), collapse=" | ")
##     }]
##     call[, mystery := !missedj & !complex & mate.mappable & seed.mappable & !insertion]
##     if(any(recall)){
##         call[recall]$call = update.call(call[recall])
##         call[(missedj), seed.mappable := TRUE]
##         call[(missedj), mate.mappable := TRUE]
##         filtered.contigs = grbind(filtered.contigs, recall.res$contigs, fill = TRUE, use.names = TRUE)
##     }
##     ## if(call$mystery){
##     ##     if(verbose) message("mystery: repeating assembly at larger intervals")
##     ##     pp = gr.stripstrand((gr.tile(wholseed-250, 500)+250) %Q% (width == 1e3))
##     ##     ## changed from ri to ri.input for consistency
##     ##     wide.contigs = .build.tigs(ri.input, pp, li$sample, li$leix, verbose = TRUE)
##     ##     if (mix.tn) {
##     ##         call = caller(li, wide.contigs, ref_obj=ref_obj, pad = 5e3)
##     ##     } else {
##     ##         call = caller(li, wide.contigs[grepl('sample', track),], ref_obj = ref_obj, pad = 5e3)
##     ##     }
##     ## }

##     call[, category := ifelse(complex, "complex rearrangement",
##                        ifelse(missedj, "missed junction",
##                        ifelse(mystery | grepl("mystery", mate.repeats) | (seed.mappable & mate.mappable),
##                               "mystery",
##                               paste("type", as.integer(!seed.mappable) + as.integer(!mate.mappable),
##                                     "loose end"))))]

##     ## update call to set category used in paper (type 0 loose ends)
##     call[, paper_category := ifelse(category == "complex rearrangement" | category == "missed junction",
##                                     "type 0 loose end",
##                                     category)]
    
##     call[category=="mystery", mystery := TRUE]
    
##     call[,  ":="(
##         leix = li[, leix],
##         loose.end = li[, paste0(seqnames, ":", start, strand)],
##         sample = id)]

##     if (verbose) {
##         message("category: ", call$category)
##     }

##     res = list(all.contigs = all.contigs,
##                call = call,
##                filtered.contigs = as.data.table(filtered.contigs))
##     ## if (exists('wide.contigs')) {
##     ##     if (inherits(wide.contigs, "data.table")) {
##     ##         res$wide.contigs = wide.contigs
##     ##     }
##     ## }
##     return(res)
## }
