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

#' @name call_loose_end_wrapper
#' @title call_loose_end_wrapper
#'
#' @description
#'
#' call loose end given reads
#'
#' @param le.dt (from prep_loose_ends)
#' @param reads.dt (from loosereads)
#' @param concat.bwa (BWA from RSeqLib) concatenated fasta
#' @param human.bwa (BWA from RSeqLib)
#' @param pad (window for local assembly)
#' @param mix.tn (logical) mix tumor and normal reads? default TRUE
#' @param minimap (logical) default FALSE
#' @param max.reads (numeric) max reads per loose end before downsampling (default 5000)
#' @param verbose
#'
#' @return list with three items:
#' $call table with calls
#' $filtered.contigs data.table with filtered contigs
#' $all.contigs data.table with all contigs
#'
#' @export
call_loose_end_wrapper = function(id = "",
                                  le.dt = data.table(),
                                  reads.dt = data.table(),
                                  concat.bwa = NULL,
                                  human.bwa = NULL,
                                  concat.fn = NA_character_,
                                  pad = 5e3,
                                  mix.tn = TRUE,
                                  minimap = FALSE,
                                  outdir = "./",
                                  max.reads = 5000,
                                  verbose = FALSE) {

    if (is.null(concat.bwa)) { stop("Must supply BWA object as concat.bwa") }
    if (is.null(human.bwa)) { stop("Must supply BWA object as human.bwa") }
        
    empty.res = list(call = data.table(), filtered.contigs = data.table(), all.contigs = data.table())

    le.dt = prep_loose_ends(li = le.dt, id = id)

    if (!le.dt[, .N]) {
        return(empty.res)
    }

    if (!reads.dt[, .N]) {
        return(empty.res)
    }

    full.res = lapply(1:le.dt[, .N],
                      function(ix) {
                          sel = dt2gr(reads.dt) %^% (dt2gr(le.dt[ix,]) + pad)
                          if (any(sel)) {
                              qns = reads.dt[sel, qname]
                              this.reads.dt = reads.dt[qname %in% qns,]
                          } else {
                              this.reads.dt = reads.dt
                          }
                          ri = prep_loose_reads(li = le.dt[ix,],
                                                loose.reads.dt = this.reads.dt)

                          sub.res = call_loose_end(li = le.dt[ix,],
                                                   ri = ri,
                                                   concat.bwa = concat.bwa,
                                                   human.bwa = human.bwa,
                                                   pad = pad,
                                                   mix.tn = mix.tn,
                                                   concat.fn = concat.fn,
                                                   minimap = minimap,
                                                   outdir = paste0(outdir, "/", ix),
                                                   max.reads = max.reads,
                                                   verbose = verbose)
                          return(sub.res)
                      })

    full.call = lapply(full.res, function(x) {x$call}) %>%
        rbindlist(fill = TRUE, use.names = TRUE)
    full.filtered = lapply(full.res, function(x) {x$filtered.contigs}) %>%
        rbindlist(fill = TRUE, use.names = TRUE)
    full.contigs = lapply(full.res, function(x) {x$all.contigs}) %>%
        rbindlist(fill = TRUE, use.names = TRUE)

    return(list(call = full.call, filtered.contigs = full.filtered, all.contigs = full.contigs))
}

#' @name build.from.win
#'
#' assembles strand-specific reads and their mates into contigs and counts read support per contig within a single seed window
#' @param win GRanges seed window
#' @param ri data.table reads and their mates
#' @param tracks optional, character which tracks (sample + strand) to assemble, default all
#' @param align.thres (numeric) width of alignment / number of non-N characters
build.from.win = function(win, ri, tracks=NULL, verbose = FALSE, align.thres = 0.9){
    if(!("track" %in% colnames(ri))){
        ri[, track := ifelse(strand=="+", "for", "rev")]
    }
    if(is.null(tracks)) tracks = ri[, unique(track)]

    if (verbose) {
        message("Tracks used for assembly: ")
        message(paste(tracks, collapse = '\n'))
    }
    ri[track %in% tracks, {
        t = track[1]
        ri[, seed := dt2gr(ri) %N% win > 0 & track==t]
        rtmp = ri[qname %in% qname[seed]]
        rtmp[(seed)][is.dup(qname), seed := 1:.N == 1, by=qname] ## in case of a foldback where both reads are "seed", choose one as the anchor
        rtmp[, seed.frame := ifelse(seed, reading.frame, as.character(reverseComplement(DNAStringSet(reading.frame))))]
        if(nrow(rtmp) < 5){
            data.table(peak = gr.string(win[,c()]), seq = as.character(NA), good.assembly=FALSE, cov=ri[(seed), .N], mapq60=ri[(seed) & mapq==60, .N], unassembled.reads=nrow(rtmp))
        } else{
            srf = rtmp[, seed.frame]
            
            if (verbose) {
                message("Starting assembly with ", length(srf), " reads")
            }
            f = Fermi(srf, assemble=T)
            contigsf = contigs(f)
            unassembled.reads = length(srf)

            if (verbose) {
                message("Initial number of contigs: ", length(contigsf))
            }
            
            if (length(contigsf)) {
                reassemble = TRUE
                while(reassemble){

                    ## using current contigs, count the number of unassembled reads
                    ctigs = BWA(seq=contigsf)
                    aln = ctigs[srf]

                    ## qnames not included in alignment
                    s1 = which(!(seq_along(srf) %in% as.integer(aln$qname)))

                    ## qnames with "bad" alignments (e.g. alignment width < 0.9 * #non-N characters)
                    aln.dt = as.data.table(aln)
                    align.frac.dt = aln.dt[, .(best.aln = max(width / nchar(gsub("N", "", seq)))), by = qname]
                    s2 = align.frac.dt[best.aln < align.thres, as.integer(qname)]
                    ## s2 = as.data.table(aln)[, max(width / nchar(gsub("N", "", seq))), by=qname][V1 < 0.9, as.integer(qname)]
                    ## if the number of unassembled reads does not change from the previous iteration...
                    if(length(s1) + length(s2) == unassembled.reads | length(s1) + length(s2) < 7) {
                        reassemble = FALSE
                    } else {
                        unassembled.reads = length(s1) + length(s2)

                        if (verbose) {
                            message("Reassembling. Number of unassembled reads: ", unassembled.reads)
                        }
                        ## if(unassembled.reads > 7){
                        ## message("Starting assembly with ", unassembled.reads, " reads")
                        c2 = Fermi(srf[c(s1, s2)], assemble=T)
                        if(length(contigs(c2))){
                            contigsf = c(contigsf, contigs(c2))
                        } else {
                            reassemble = FALSE
                        }
                    }
                }
                goodassembly = rep(TRUE, length(contigsf))
                for(i in seq_along(contigsf)){
                    ctig = BWA(seq = contigsf[i])
                    ra = ctig[srf]
                    sc = table(factor(strand(ra), c("+", "-")))
                    if(!any(sc[1:2] == 0)) goodassembly[i] = FALSE
                    if(sc["+"] < sc["-"]) contigsf[i] = as.character(reverseComplement(DNAStringSet(contigsf[i])))
                }
                data.table(peak = gr.string(win[,c()]), seq = contigsf, good.assembly=goodassembly, cov=ri[(seed), .N], mapq60=ri[(seed) & mapq==60, .N], unassembled.reads=unassembled.reads)
            } else{
                data.table(peak = gr.string(win[,c()]), seq = as.character(NA), good.assembly=FALSE, cov=ri[(seed), .N], mapq60=ri[(seed) & mapq==60, .N], unassembled.reads = length(srf))
            }
        }
    }, by=track]
}

#' match.seq
#'
#' returns logical vector of length subject indicating whether matches were found with any query
#' @param query PDict of query sequences to match
#' @param subject DNAStringSet of sequences to parse for matches

match.seq = function(query, subject)
{
    if (is.null(names(subject)))
        names(subject) = 1:length(subject)

    if (any(duplicated(names(subject))))
    {
        warning('Names of subject sequences have duplicates, deduping')
        names(subject) = dedup(names(subject))
    }

    totals = base::lengths(Biostrings::vwhichPDict(query, subject)) > 0
    return(totals)
}

#' eighteenmer
#'
#' generate a PDict of telomeric 18mers
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

#' munch
#' 
#' returns logical vector of length reads
#' indicating does read contain any match to query
#' searches for matches on same strand as input seq
#' @param reads GRanges or data.table containing $seq field with character sequence
#' @param query optional, PDict of motifs, default= 18mers of telomeric motifs (both strands)
#' @export
munch = function(reads, query=NULL) {
    if(is.null(query)) query = eighteenmer()
    seq = reads$seq
    seq_frame_ref = DNAStringSet(seq)
    return(match.seq(query, seq_frame_ref))
}

#' @name return.dt
#' @title return.dt
#'
#' @description
#'
#' construct loose end call from a set of binary inputs
#'
#' @param reference (logical) reference alignment?
#' @param complex (logical) complex rearrangement?
#' @param missedj (logical) type 0 loose end
#' @param novel (logical) contains unaligned sequence?
#' @param mystery (logical) either seed or mate is mystery
#' @param insertion (logical) contig contain insertion?
#' @param rrep (logical) type of reference (human) repeat
#' @param irep (logical) non-reference (non-human) repeat
#' @param refmap (logical) seed alignment is mappable
#' @param novmap (logical) mate alignment is mappable
#' @param seed.only (logical) contig overlaps seed region completely
#' @param mate.only (logical) contig does not overlap seed region at all
#' @param junction (character)
#' constructs output data.table from logical arguments
return.dt = function(reference, complex, missedj, novel, mystery, insertion, refmap=NULL, novmap=NULL, rrep=NULL, irep=NULL, nrep=NULL, seed.only = NA, mate.only = NA, junction=""){
    if(length(rrep)) if(length(rrep %Q% (c_spec == "unassembled"))) rrep[rrep$c_spec=="unassembled"]$c_spec = "low MAPQ"
    if(length(nrep)) if(length(nrep %Q% (c_spec == "unassembled"))) nrep[nrep$c_spec=="unassembled"]$c_spec = "low MAPQ"
    if(length(irep)) if(length(irep %Q% (c_spec == "unassembled"))) irep[irep$c_spec=="unassembled"]$c_spec = "low MAPQ"
    if(length(rrep)) if(length(rrep %Q% (c_spec != "low MAPQ"))) rrep = rrep %Q% (c_spec != "low MAPQ")
    if(length(nrep)) if(length(nrep %Q% (c_spec != "low MAPQ"))) nrep = nrep %Q% (c_spec != "low MAPQ")
    if(length(irep)) if(length(irep %Q% (c_spec != "low MAPQ"))) irep = irep %Q% (c_spec != "low MAPQ")

    frep = c(irep, nrep)
    rreps = ifelse(length(rrep), (rrep %Q% (rev(order(width))))[1]$c_spec, "")
    ireps = ifelse(length(irep), (irep %Q% (rev(order(width))))[1]$c_spec, "")
    nreps = ifelse(length(nrep), (nrep %Q% (rev(order(width))))[1]$c_spec, "")
    freps = ifelse(length(frep), (frep %Q% (rev(order(width))))[1]$c_spec, "")
    if(mystery){
        nreps = freps = "mystery"
        novmap=FALSE
    }
    call.string = NULL
    if(missedj) call.string = "missed junction"
    if(complex) call.string = paste(c(call.string, "complex"), collapse="; ")
    if(reference){
        radd = paste(ifelse(refmap, "mappable", "unmappable"), "seed repeat")
        if(rreps == "") rreps = "low MAPQ"
        radd = paste0(radd, ":", rreps)
        call.string = paste(c(call.string, radd), collapse="; ")
    }
    if(novel){
        nadd = paste(ifelse(novmap, "mappable", "unmappable"), "mate repeat")
        if(nreps == "") nreps = "low MAPQ"
        nadd = paste0(nadd, ":", nreps)
        call.string = paste(c(call.string, nadd), collapse="; ")
    }
    if(insertion){
        iadd = "insertion"
        if(ireps != "") iadd = paste0(iadd, ":", ireps)
        call.string = paste(c(call.string, iadd), collapse="; ")
    }
    if(mystery) call.string = paste(c(call.string, "mystery"), collapse="; ")
    if(is.null(refmap)) refmap = logical()
    if(is.null(novmap)) novmap = logical()
    ##    single.call = ifelse(complex, "complex rearrangement", ifelse(missedj, "missed junction", ifelse(reference & !(refmap), "seed repeat", ifelse(novel & !(novmap), "mate repeat", "mystery"))))
    single.call = ifelse(complex, "complex rearrangement", ifelse(missedj, "missed junction", ifelse(mystery | grepl("mystery", freps), "mystery", paste("type", as.integer(!refmap) + as.integer(!novmap), "loose end"))))
    if(single.call == "mystery") mystery = TRUE
    data.table(seedrep=reference,
               complex=complex,
               missedj=missedj,
               materep=novel,
               mystery=mystery,
               insertion=insertion,
               seed.repeats = rreps,
               ins.repeats = ireps,
               nov.repeats = nreps,
               mate.repeats = freps,
               seed.mappable = refmap,
               mate.mappable = novmap,
               junction = junction,
               call = call.string,
               category = single.call,
               seed.only = seed.only,
               mate.only = mate.only)
}

#' update.call
#'
#' modifies call string if logical column values have changed
#' (used when supplementing contig calls with discordant read pairs)
update.call = function(dt){
    dt = copy(dt)
    dt[, i := 1:.N]
    call.string = dt[, {
        call.string = character()
        if(missedj) call.string = "missed junction"
        if(complex) call.string = paste(c(call.string, "complex"), collapse="; ")
        if(seedrep){
            radd = paste(ifelse(seed.mappable, "mappable", "unmappable"), "seed repeat")
            radd = paste0(radd, ":", seed.repeats)
            call.string = paste(c(call.string, radd), collapse="; ")
        }
        if(materep){
            nadd = paste(ifelse(mate.mappable, "mappable", "unmappable"), "mate repeat")
            nadd = paste0(nadd, ":", nov.repeats)
            call.string = paste(c(call.string, nadd), collapse="; ")
        }
        if(insertion){
            iadd = "insertion"
            if(ins.repeats != "") iadd = paste0(iadd, ":", ins.repeats)
            call.string = paste(c(call.string, iadd), collapse="; ")
        }
        if(mystery) call.string = paste(c(call.string, "mystery"), collapse="; ")
        call.string
    }, by=i]
    return(call.string$V1)
}

#' merge.telomeres
#'
#' construct a GRanges in contig coordinates for telomere motif matches
merge.telomeres = function(treps, x, ch, out.dt, c_spec = "telomere"){
    values(treps) = out.dt[.(as.character(seqnames(treps))), colnames(ch@values), with=F]
    treps$AS = width(treps)
    treps$y = gr.string(dt2gr(as.data.table(treps)[, seqnames := "tel"]))
    treps$cigar = ifelse(start(treps) > 1, paste0(start(treps)-1, "S", width(treps), "M"), paste0(width(treps), "M"))
    treps$c_type = "rep"
    treps$mapq = 0
    treps$c_spec = c_spec
    treps$seed = FALSE
    x = dt2gr(rbind(as.data.table(x), as.data.table(treps), fill=TRUE, use.names=TRUE))
    return(x)
}

#' find.telomeres
#'
#' search for telomeric motif matches in contigs
find.telomeres = function(query, out.dt){
    gd = munch(out.dt, query)
    if(!(any(gd))) return(NULL)
    qns = out.dt[gd, qname]
    do.call('c', lapply(qns, function(qn) GenomicRanges::reduce(GRanges(qn, IRanges(do.call('c', lapply(Biostrings::vwhichPDict(query, DNAStringSet(out.dt[.(qn), seq]))[[1]], function(i) as.integer(gregexpr(as.character(query@dict0)[i], out.dt[.(qn), seq])[[1]]))), width=18), strand="+"))))
}

#' caller
#'
#' parses contig alignments to assign loose end to category and describe repeat types
#' @param li data.table loose end to evaluate
#' @param calns optional, data.table of contig alignments, default=NULL (will not parse alignments)
#' @param insert optional, integer pad representing insert length in bp to identify seed alignments, default=750
#' @param pad optional, window around loose end to allow contig seed windows, default=1e3
#' @param uannot optional, GRanges of unmappable annotations, default bin/101.unmappable.annotations.rds
#' @param ref_dir optional, path to directory of unzipped reference tarball, default assumes 'package/extdata/hg19_looseends'
#' @param ref_obj optional, list of BWA objects built from ref_dir fastas, names must match expected "human" "rep" "polyA" "microbe" (only "human" is used), default=NULL
#' @param return.contigs (logical)
#' @export
caller = function(li,
                  calns = NULL,
                  insert = 750,
                  pad=NULL,
                  uannot=NULL,
                  human = NULL,
                  return.contigs = FALSE) {
    reference = FALSE
    complex = FALSE
    missedj = FALSE
    novel = FALSE
    mystery = FALSE
    insertion = FALSE
    refmap = FALSE
    novmap = FALSE
    rrep = nrep = irep = GRanges(c_spec=character())
    junction = ""
    if(is.null(pad)) pad = 1e3
    if(is.null(uannot)){
        uannot = readRDS(system.file('extdata', '101.unmappable.annotations.rds', package='loosends'))
    }
    ## if(!is.null(ref_obj)) { human = ref_obj$human
    ## } else{
    ##     if(!file.exists(ref_dir)) stop("Provide correct ref_dir containing reference .fa files")
    ##     if(!file.exists(paste(ref_dir, "human_g1k_v37_decoy.fasta", sep="/"))) stop("ref_dir must contain human_g1k_v37_decoy.fasta")
    ##     human = BWA(fasta=paste(ref_dir, "human_g1k_v37_decoy.fasta", sep="/"))
    ## }

    if (is.null(human)) { stop("Human cannot be null, please supply BWA object") }
    
    if(is.null(calns) | !nrow(calns)) {
        nrep = uannot[1]
        nrep$c_spec = "mystery"
        if(dt2gr(li) %N% uannot){
            reference = TRUE
            rrep = uannot %&% dt2gr(li)
            if(is.na(rrep$c_spec)) rrep$c_spec = rrep$repClass
        } else {
            refmap = TRUE
        }
        mystery = TRUE
        if (return.contigs) {
            res = list(call = return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, irep=irep, nrep=nrep, refmap=refmap, novmap=novmap),
                       irep = irep,
                       nrep = nrep,
                       contigs = GRanges())
            return(res)
        }
                       
        return(return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, irep=irep, nrep=nrep, refmap=refmap, novmap=novmap))
    }
    if(inherits(calns$mapq, "character")) calns$mapq = as.integer(calns$mapq)
    calns[c_spec == "polyA", c_type := "rep"]
    calns[c_spec == "ribosomal", c_type := "ribosomal"]
    seed = GRanges(calns$peak)
    strand(seed) = ifelse(grepl("for", calns$track), "+", "-")
    seed$qname = calns$qname
    ## only use contigs overlapping loose end within 1kb window 
    if("seed" %in% colnames(calns)) calns$seed = NULL
    ## we should ignore loose end strand here?
    calns = calns[!is.na(gr.match(seed, gr.flipstrand(dt2gr(li))+pad, ignore.strand=T))]
    if(!nrow(calns)){
        nrep = uannot[1]
        nrep$c_spec = "mystery"
        if(dt2gr(li) %N% uannot){
            reference = TRUE
            rrep = uannot %&% dt2gr(li)
            if(is.na(rrep$c_spec)) rrep$c_spec = rrep$repClass
        } else {
            refmap = TRUE
        }
        mystery = TRUE
        if (return.contigs) {
            res = list(call = return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, irep=irep, nrep=nrep, refmap=refmap, novmap=novmap),
                       irep = irep,
                       nrep = nrep,
                       contigs = GRanges())
            return(res)
        }
        return(return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, irep=irep, nrep=nrep, refmap=refmap, novmap=novmap))
    }
    
    ## filter out contigs with too few good input reads?
    ## for now, we won't use MAPQ as a filter
    if(calns[map60.cov >= 0, .N == 0]){
        reference = TRUE
        rrep = GRanges(li$seqnames, IRanges(li$start, li$end), c_spec="low MAPQ")
    }

    ## checking this isn't just the seed region fully mapping & not going anywhere
    ## also making sure it doesn't leave but come back
    seed = GRanges(calns$peak)
    strand(seed) = ifelse(grepl("for", calns$track), "+", "-")
    seed$qname = calns$qname
    ## if within insert (default 750) of the seed of the contig, on the same strand
    ## then mark as a seed alignment
    calns$seed = !is.na(gr.match(dt2gr(calns), seed+insert, by="qname", ignore.strand=T))
    if(any(calns$seed) & all(calns[(seed), mapq == 60])) refmap = TRUE
    ch = cgChain(calns)
    x = ch@.galx
    values(x) = ch@values
    x$y = gr.string(dt2gr(calns)[floor(as.numeric(rownames(ch@values)))])
    ## ignoring indel cigar alignments -- w/r/t when width is used to prioritize repeats
    x$source.id = floor(as.numeric(rownames(ch@values)))
    ctl = si2gr(x)
    x = dt2gr(as.data.table(x)[, ":="(start = min(start), end = max(end)), by=source.id][!duplicated(source.id)][,!"source.id"])
    seqlengths(x) = seqlengths(ctl)[seqlevels(x)]
    
    ## now that we're in contig coordinates, parse for telomeres
    out.dt = calns[!duplicated(qname)]
    setkey(out.dt, qname)
    tgreps = find.telomeres(eighteenmer('g'), out.dt)
    tcreps = find.telomeres(eighteenmer('c'), out.dt)
    if(!is.null(tgreps)) x = merge.telomeres(tgreps, x, ch, out.dt, c_spec="G telomere")
    if(!is.null(tcreps)) x = merge.telomeres(tcreps, x, ch, out.dt, c_spec="C telomere")
    seqlengths(x) = seqlengths(ctl)[seqlevels(x)]
    ## viral calls overlapping telomere matches are ignored -- viral reference seqs include (TTAGGG)n
    ## actually viral calls overlapping anything human are ignored??
    if(any(x$c_type=="viral") & (!is.null(tgreps) | !is.null(tcreps))){
        virs = x %Q% (c_type == "viral")
        kvir = virs %O% (x %Q% (grepl("telomere", c_spec)))
        x = c(x %Q% (c_type != "viral"), virs[kvir < 0.9])
    }
    ## ribosomal calls overlapping polyA calls are ignored -- ribosomal ref includes (A)n
    if(any(x$c_spec=="ribosomal") & any((x %Q% (c_spec=="ribosomal")) %N% (x %Q% (c_spec=="polyA")))){
        ribs = x %Q% (c_spec == "ribosomal")
        krib = ribs %O% (x %Q% (c_spec == "polyA"))
        x = c(x %Q% (c_spec != "ribosomal"), ribs[krib < 0.9])
    }
    ## prioritze telomere motif matches over telomeric coordinate alignments
    if(any(x$c_spec=="Telomere / sub") & (!is.null(tgreps) | !is.null(tcreps))){
        tela = x %Q% (c_spec == "Telomere / sub")
        utel = is.na(gr.match(tela, x %Q% (grepl("tel", y)), by="seed"))
        x = c(x %Q% (c_spec != "Telomere / sub"), tela[utel])
    }

    
    ## make sure the contig goes and stays away:
    ## last 20 bases of the contig cannot align to seed window
    ends = gr.end(ctl, width=20)
    ends$seed = ends %N% (x %Q% (seed))
    ends$human60 = ends %N% (x %Q% (c_spec=="human") %Q% (mapq==60) %Q% (!(seed)))
    ## outp: qnames of contigs with a seed alignment in the last 20 bases
    outp = as.data.table(ends)[seed > 0, as.character(seqnames)]

    if(reference){
        rrep = copy(x[1])
        rrep$c_spec = "low MAPQ"
    }

    ## add to reference repeats: contig seed regions that also align elsewhere, at least 20bp
    ## (should be making sure the overlap specifically is at least 20bp)
    rrep = c(rrep, x %&%
                   (x %Q% (seed)) %Q% (mapq<60 | c_type!="human") %Q% (width > 19) %Q% (seqnames %in% outp))
    if(any(rrep$c_spec == "human" & rrep$mapq < 60))
        rrep[rrep$c_spec == "human" & rrep$mapq < 60]$c_spec = "low MAPQ"
    rrep = rrep %Q% (c_spec!="human")
    if(length(rrep)) {
        reference = TRUE
    }

    ## filter out contigs with seed alignment in the last 20 bases
    x2 = x %Q% (!(seqnames %in% outp)) ## save x so you can check for reference repeats
    if(!length(x2)){
        nrep = uannot[1]
        nrep$c_spec = "mystery"
        if(dt2gr(li) %N% uannot){
            reference = TRUE
            sup = uannot %&% dt2gr(li)
            if(is.na(sup$c_spec)) sup$c_spec = sup$repClass
            rrep = grbind(rrep, sup)
            rrep$c_spec = as.character(rrep$c_spec)
        } else refmap = refmap | all((x %Q% (seed))$mapq == 60)
        mystery = TRUE
        if (return.contigs) {
            res = list(call = return.dt(reference, complex,
                                        missedj, novel, mystery,
                                        insertion, rrep=rrep,
                                        irep=irep,
                                        nrep=nrep,
                                        refmap=refmap,
                                        novmap=novmap,
                                        seed.only = TRUE),
                       irep = irep,
                       nrep = nrep,
                       contigs = x2)
            return(res)
        }
        return(return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, irep=irep, nrep=nrep, refmap=refmap, novmap=novmap))
    } else x = x2
    
    ## trim the lowMQ bit to be only outside the repeat
##    if(any(x$c_type=="rep") & any((x %Q% (c_type=="rep")) %N% (x %Q% (c_spec=="human") %Q% (!seed & mapq < 60)))){
##        mq = x %Q% (!seed & c_spec == "human") %Q% (mapq < 60)
##        nmq = GenomicRanges::setdiff(mq, x %Q% (c_type=="rep"))
##        ov = gr.findoverlaps(nmq, mq, ignore.strand=F)
##        values(ov) = values(mq[ov$subject.id])
##        x = c(x %Q% (seed | c_spec != "human" | mapq==60), ov)
##    }

    ## remember you won't always have caught the junction -- might just be the other side
    sections = disjoin(x)
    sections$score = sections %N% x
    sections$seeds = sections %N% (x %Q% (seed))
    sections$seed60 = sections %N% (x %Q% (seed & mapq==60))
    sections$hg19 = sections %N% (x %Q% (c_type=="human"))
    sections$human = sections %N% (x %Q% (c_spec=="human"))
    sections$human60 = sections %N% (x %Q% (c_spec=="human") %Q% (mapq==60))

    x = x %**% sections
    if(any(x$c_type=="viral")){
        virs = x %Q% (c_type == "viral")
        vwid = as.data.table(virs)[, sum(width), by=.(query.id, seed, qname)]
        nh = (x %Q% (c_type != "human" | mapq < 60)) %Q% (c_type != "viral")
        nwid = as.data.table(nh)[, sum(width), by=.(query.id, seed,qname)][, max(V1), by=.(seed, qname)]
        nwid[, sq := paste(seed, qname)]
        setkey(nwid, sq)
        vwid[, nwid := nwid[.(paste(vwid$seed, vwid$qname)), V1]]
        vwid[is.na(nwid), nwid := 0]
        virs = virs %Q% (query.id %in% vwid[V1 > nwid][, query.id])
        x = c(x %Q% (c_type != "viral"), virs)

        sections = disjoin(x)
        sections$score = sections %N% x
        sections$seeds = sections %N% (x %Q% (seed))
        sections$seed60 = sections %N% (x %Q% (seed & mapq==60))
        sections$hg19 = sections %N% (x %Q% (c_type=="human"))
        sections$human = sections %N% (x %Q% (c_spec=="human"))
        sections$human60 = sections %N% (x %Q% (c_spec=="human") %Q% (mapq==60))
    }
    
    x = as.data.table(x)
    x[, clength := seqlengths(ctl)[qname]]

    if(x[, any(seed)] & x[(seed), all(mapq<60)]){
        reference = TRUE
        rrep = rrep %Q% (!(seqnames %in% outp))
        if(length(setdiff(colnames(values(dt2gr(x))), colnames(values(rrep)))))
            rrep = rrep %$% dt2gr(x)[, setdiff(colnames(values(dt2gr(x))), colnames(values(rrep)))]
        rrep = grbind(rrep, dt2gr(copy(x[(seed)][!duplicated(qname)])[, c_spec := "low MAPQ"], seqlengths=seqlengths(ctl)))
        rrep$c_spec = as.character(rrep$c_spec)
    }
    x = dt2gr(x)
    seqlengths(x) = seqlengths(ctl)[seqlevels(x)]
    sections = sections %$% (x[, 'clength']) %$% (x %Q% (c_spec != "human"))[, 'c_spec']
    sections = as.data.table(sections)
    
    sections[, seed.end := ifelse(any(seeds>0), max(end[seeds>0], na.rm=T), as.integer(0)), by=seqnames]
    sections[, seed.side := end <= seed.end]
    sections[, mate.side := start > seed.end]

    ## require at least 20 bp of seed and mate
    ## in addition, at least 20 bp need to align with MAPQ60 to be considered mappable
    sections[, seed.map := ifelse(seed.end<19,
                                  "mystery",
                           ifelse(##any(seed60 > 0 & seed.side & width > 19, na.rm = TRUE),
                                  any(seed60 == seeds & seed60 > 0 & seed.side & width>19),
                                  "mappable",
                                  "unmappable")), by=seqnames]
    sections[is.na(seed.map), seed.map := "mystery"]
    sections[, mate.map := ifelse(seed.end+19 >= clength, "mystery",
                           ifelse(## any(human60 > 0 & mate.side & width > 19, na.rm = TRUE),
                                  sum(width[human60 > 0 & mate.side]) >= clength-seed.end-19,
                                  "mappable",
                                  "unmappable")), by=seqnames]
    sections[is.na(mate.map), mate.map := "mystery"]
    if(length(ends %Q% (seed == 0) %Q% (human60 > 0))){
        me = as.character(seqnames(ends %Q% (seed == 0) %Q% (human60 > 0)))
        sections[seqnames %in% me, mate.map := ifelse(sum(width[human60 > 0 & mate.side]) >= clength-seed.end-99, "mappable", "unmappable"), by=seqnames]
    }

    ## keep contigs where seeds and mates both are non-mystery
    sm = sections[, any(seed.map!="mystery" & mate.map!="mystery"), by=seqnames][(V1), as.character(seqnames)]
    if(length(sm)){
        sections = sections[seqnames %in% sm]
        x = x %Q% (seqnames %in% sm)
        seqlevels(x) = sm
    }

    ## if > 99 bp of the contig is on the seed side but does not overlap the peak
    ## then classify as mystery...
    ## since this implies that there is a big insertion in the seed region
    cs = sections[, sum(width[seed.side & !seeds]) > 99, by=seqnames][(V1), as.character(seqnames)]
    if(length(cs)){
        if(all(sections[, seqnames %in% cs])){
            mystery = TRUE
            nrep = dt2gr(sections[(mate.side)][1][, c_spec := "mystery"])
            novmap = FALSE
            if (return.contigs) {
                res = list(call = return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, irep=irep, nrep=nrep, refmap=refmap, novmap=novmap),
                           irep = irep,
                           nrep = nrep,
                           contigs = x)
                return(res)
            }
            return(return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, irep=irep, nrep=nrep, refmap=refmap, novmap=novmap, junction=junction))
        }
        sections = sections[!(seqnames %in% cs)]
        x = x %Q% (!(seqnames %in% cs))
        seqlevels(x) = setdiff(seqlevels(x), cs)
    }

    ## keep track of if contig is seed only or mate only
    mate.only = sections[,all(seed.map=="mystery")]
    seed.only = sections[,all(mate.map=="mystery")]
    
    if(sections[,all(seed.map=="mystery")]){
        if(dt2gr(li) %N% uannot){
            reference = TRUE
            sup = (uannot %&% dt2gr(li) %Q% (width == max(width)))[1]
            rrep = x %Q% (!duplicated(qname))
            rrep$c_spec = sup$c_spec
            rrep$c_type = "rep"
            if(is.na(sup$c_spec)) rrep$c_spec = sup$repClass
            sections[, seed.map := "unmappable"]
        } else{
            cs = out.dt[qname %in% as.character(sections$seqnames), GRanges(paste0("chr", peak), seqinfo=seqinfo(Hsapiens), qname=qname)]
            cs$seq = as.character(getSeq(Hsapiens, cs))
            ca = human[cs$seq] %Q% (rev(order(mapq))) %Q% (!duplicated(qname))
            cm = setNames(as.integer(ca$mapq), cs[as.integer(ca$qname)]$qname)
            sections[, seed.map := ifelse(cm[as.character(seqnames)]==60, "mappable", "unmappable")]
            if(any(cm < 60)){
                sup = out.dt[qname %in% names(cm)[cm<60]]
                sup$c_type = "human"
                sup$c_spec = "low MAPQ"
                rrep = c(rrep, dt2gr(sup, seqlengths = seqlengths(x)))
            }
        }
    }    

    sections = dt2gr(as.data.table(sections))
    seqlengths(sections) = seqlengths(ctl)[seqlevels(sections)]
    

    ## can only catch a reference repeat this way if the seed is part of the contig
    rr = unique(as.character(seqnames(x %Q% (seeds>0) %Q% (score > 1 | human>human60) %Q% (width > 19))))
    if(length(rr)){
        reference = TRUE
        rrep = rrep %Q% (seqnames %in% as.character(seqnames(x)))
        if(length(setdiff(colnames(values(x)), colnames(values(rrep)))))
            rrep = rrep %$% x[, setdiff(colnames(values(x)), colnames(values(rrep)))]
        rrep = c(rrep, x %Q% (seqnames %in% rr) %Q% (seeds>0) %Q% (!(seed) | seeds > 1 | human>human60) %Q% (width > 19))
        rrep = rrep %Q% (!(c_spec=="human" & mapq==60))
        if(any(rrep$c_spec == "human" & rrep$mapq < 60))
            rrep[rrep$c_spec == "human" & rrep$mapq < 60]$c_spec = "low MAPQ"
    }    

    ## make sure if there are multiple they go to the same place..
    ## seed has good MAPQ -- if not included in any contig, align using getSeq
    ## last 20 bases overlap good MAPQ
    ## no MAPQ60 gaps 100bp+
    ## first check gaps:
    gapcs = gaps(sections) %Q% (strand=="+") %Q% (width > 19)
    sections$gap = as.character(seqnames(sections)) %in% as.character(seqnames(gapcs))
    if(any(sections$gap & sections$mate.map=="mappable")) sections[sections$gap & sections$mate.map=="mappable"]$mate.map = "unmappable"
    
    mj = setdiff(levels(seqnames(x)), unique(as.character(seqnames(gaps(sections %Q% (human60>0)) %Q% (strand == "+") %Q% (width > 99)))))
    mj = as.data.table(sections)[seqnames %in% mj, any(seed.map=="mappable" & mate.map=="mappable"), by=seqnames][(V1), as.character(seqnames)]
    ## next check end good MAPQ:
    if(length(mj))
        mj = as.character(seqnames(ends %Q% (seqnames %in% mj) %Q% (human60 > 0)))
    ## finally check seed qual:
    if(length(mj)){
        xs = x %Q% (seed) %Q% (rev(order(mapq))) %Q% (!duplicated(seqnames))
        cm = setNames(xs$mapq, xs$qname)
        if(length(cm) & all(cm < 60)) mj = mj[0]
        if(any(!(mj %in% names(cm)))){
            cs = GRanges(paste0("chr", out.dt[qname %in% mj][!(qname %in% names(cm))]$peak), seqinfo=seqinfo(Hsapiens))
            cs$qname = out.dt[qname %in% mj][!(qname %in% names(cm))]$qname
            cs$seq = as.character(getSeq(Hsapiens, cs))
            ca = human[cs$seq] %Q% (rev(order(mapq))) %Q% (!duplicated(qname))
            cm = c(cm, setNames(as.integer(ca$mapq), cs[as.integer(ca$qname)]$qname))
        }
        mj = mj[cm[mj] == 60]
    }
    if(length(mj)){
        missedj = TRUE
        refmap = TRUE
        novmap = TRUE

        rrep = rrep %Q% (!is.na(seqnames)) %Q% (seqnames %in% mj)
        reference = length(rrep) > 0
        nrep = nrep %Q% (!is.na(seqnames)) %Q% (seqnames %in% mj)
        novel = length(nrep) > 0
        sections = sections %Q% (seqnames %in% mj)
        seqlevels(sections) = mj
        x = x %Q% (seqnames %in% mj)
        seqlevels(x) = mj

        mate = GenomicRanges::reduce(parse.gr((x %Q% (seqnames %in% mj) %Q% (!seed) %Q% (mapq==60) %Q% (c_spec == "human"))$y) + 100)
        junction = paste(c(gr.string(dt2gr(li)), gr.string(gr.start(mate, ignore.strand=F))), collapse = " | ")
        ## can't just be one -- in case of foldbacks break glass
        if(length(mate)>1){
            if(as.data.table(x %Q% (seqnames %in% mj) %Q% (!seed) %Q% (mapq==60) %Q% (c_spec == "human"))[width > 9][,.N,by=qname][,any(N>1)]){
                complex = TRUE
            } else if(!length((x %Q% (seqnames %in% mj) %Q% (seed) %Q% (mapq==60) %Q% (c_spec == "human")))){
                complex = TRUE
            } else{
                s = GenomicRanges::reduce(gr.flipstrand(parse.gr((x %Q% (seqnames %in% mj) %Q% (seed) %Q% (mapq==60) %Q% (c_spec == "human"))$y)))
                if(length(GenomicRanges::reduce(grbind(mate, s))) > 2){
                    complex = TRUE
                }
            }
        }
        inserts = gaps(sections %Q% (seqnames %in% mj) %Q% (human60>0 | seeds>0)) %Q% (strand == "+") %Q% (width > 19)
        if(length(inserts)){
            insertion = TRUE
            irep = x %&% inserts %Q% (c_spec!="human" | mapq < 60) %Q% (!seeds)
            if(any(irep$c_spec == "human" & irep$mapq < 60))
                irep[irep$c_spec == "human" & irep$mapq < 60]$c_spec = "low MAPQ"
##            if(length(irep)){
##                wides = as.data.table(gr.reduce(irep, by="query.id"))[, query.id[width > 19]]
##                irep = irep %Q% (c_spec!="human") %Q% (c_type!="human" | hg19==score) %Q% (query.id %in% wides) %Q% (!duplicated(c_spec))
##            }
        }
    }

    rrep = c(rrep, x %Q% (seeds>0) %Q% (c_spec!="human" | mapq < 60))
    nrep = c(nrep, x %Q% (seeds==0) %Q% (c_spec!="human" | mapq < 60))
    if(any(nrep$c_spec == "human" & nrep$mapq < 60))
        nrep[nrep$c_spec == "human" & nrep$mapq < 60]$c_spec = "low MAPQ"
    if(length(nrep)){
        if(any(is.na(nrep$query.id))){
            max.qi = max(nrep$query.id, na.rm=T)
            if(!is.finite(max.qi)) max.qi = 0
            nrep[is.na(nrep$query.id)]$query.id = 1:sum(is.na(nrep$query.id)) + max.qi
        }
        si = seqinfo(nrep)
        nrep = gr.reduce(nrep, by="query.id") %Q% (width > 19)
        seqinfo(nrep) = si[seqlevels(nrep)]
        if(length(nrep)){
            novel = TRUE
        }
    }

    if(any(rrep$c_spec == "human" & rrep$mapq < 60))
        rrep[rrep$c_spec == "human" & rrep$mapq < 60]$c_spec = "low MAPQ"
    if(length(rrep)){
        if(any(is.na(rrep$query.id))){
            max.qi = max(rrep$query.id, na.rm=T)
            if(!is.finite(max.qi)) max.qi = 0
            rrep[is.na(rrep$query.id)]$query.id = 1:sum(is.na(rrep$query.id)) + max.qi
        }
        si = seqinfo(rrep)
        rrep = gr.reduce(rrep, by="query.id") %Q% (width > 19)
        seqinfo(rrep) = si[seqlevels(rrep)]
        if(length(rrep)){
            reference = TRUE
        }
    }

    nr = as.character( seqnames(sections %Q% (seed.map == "mappable"))) %>% unique
    if(length(nr)){
        refmap = TRUE
        rrep = rrep %Q% (!is.na(seqnames)) %Q% (seqnames %in% nr)
        reference = length(rrep) > 0
        nrep = nrep %Q% (!is.na(seqnames)) %Q% (seqnames %in% nr)
        novel = length(nrep) > 0

        sections = sections %Q% (seqnames %in% nr)
        seqlevels(sections) = nr
        x = x %Q% (seqnames %in% nr)
        seqlevels(x) = nr
    }

    nh = as.character(seqnames(sections %Q% (mate.map == "mappable"))) %>% unique
    if(length(nh)){
        novmap = TRUE
        rel = sections %Q% (seqnames %in% nh)
        if(all(rel$seed.map == "mystery")){
            if(dt2gr(li) %N% uannot){
                reference = TRUE
                rrep = uannot %&% dt2gr(li)
                if(is.na(rrep$c_spec)) rrep$c_spec = rrep$repClass
            } else {
                refmap = TRUE
            }
        } else{
            refmap = any(rel$seed.map=="mappable")
        }
        nrep = nrep %Q% (!is.na(seqnames)) %Q% (seqnames %in% nh)
        novel = length(nrep) > 0
        if(any(nrep$c_spec == "human" & nrep$mapq < 60))
            nrep[nrep$c_spec == "human" & nrep$mapq < 60]$c_spec = "low MAPQ"
        if(!"human60" %in% colnames(values(nrep))) nrep = nrep %$% sections        
        if(insertion & any(!(nrep$human60))){
            nrep = nrep %Q% (human60)
        }

        sections = sections %Q% (seqnames %in% nh)
        seqlevels(sections) = nh
        x = x %Q% (seqnames %in% nh)
        seqlevels(x) = nh
    }

    g = gaps(sections) %Q% (strand == "+") %Q% (width > 19)
    if(length(g)){
        g = g[1]
        una = x[1]
        una$c_spec = "unaligned sequence"
        nrep = c(nrep, una)
        novel = TRUE
    }

    if(missedj) refmap = novmap = TRUE
    mystery = mystery | !any(missedj, complex, reference, novel)

    call.dt = return.dt(reference, complex, missedj, novel, mystery, insertion, rrep=rrep, irep=irep, nrep=nrep, refmap=refmap, novmap=novmap, junction=junction, seed.only = seed.only, mate.only = mate.only)

    if (return.contigs) {
        res = list(contigs = x, call = call.dt, irep = irep, nrep = nrep)
        return(res)
    }
    return(call.dt)
}


#' gr.sum.strand
#'
#' strand-specific implementing gr.sum
gr.sum.strand = function(gr){
    pos = gr.sum(gr %Q% (strand == "+"))
    strand(pos) = "+"
    neg = gr.sum(gr %Q% (strand == "-"))
    strand(neg) = "-"
    return(c(pos, neg))
}

#' read.based
#'
#' classifying loose end based on discordant read alignments by constructing pseudo-contigs based on consensus alignment patterns
#' @param li data.table loose end to classify
#' @param ri data.table read alignments
#' @param pad optional, padding around loose end li to search for discordant reads, default=1e3
#' @param human (human BWA)
#' @param return.contigs return contigs if TRUE, otherwise just return
#'
#' @return if return.contigs is TRUE, returns a list with elements $call, $filtered.contigs, $all.contigs. otherwise returns a data.table with calls.
read.based = function(li, ri, pad=NULL, human = NULL, return.contigs = FALSE, mapq.thresh = 60){
    if(is.null(pad)) pad = 1e3
    t = ifelse(li$strand == "+", "sample.rev", "sample.for")
    ri = ri[rev(order(qwidth))][!duplicated(paste(R1, qname))]
    pp = (gr.tile(gr.flipstrand(dt2gr(li))+pad, 200) %Q% (width == 200))[,c()]

    out.dt = rbindlist(lapply(1:length(pp), function(i){
        win = pp[i]
        ri[, seed := dt2gr(ri) %N% win > 0 & track==t]        
        rtmp = ri[qname %in% qname[seed]]
        rtmp[(seed)][is.dup(qname), seed := R1 == R1[1], by=qname] ## in case of a foldback where both reads are "seed", choose one as the anchor
        rtmp[is.dup(paste(qname, R1)), strand := ifelse(seed | !any(seed), strand, ifelse(strand=="-", "+", "-"))] ## in case of split alignment OF SEED READ, flip the strand of non-seed split (keep split of mate)
        ## ONLY USES MAPQ 60 READS FOR CALLING??
        rtmp = rtmp[mapq >= mapq.thresh]
        if(nrow(rtmp)==0 | rtmp[(seed), .N==0] | rtmp[!(seed), .N==0]) return(NULL)
        ctig = GenomicRanges::reduce(gr.sum.strand(dt2gr(rtmp[(seed)])) %Q% (score > 0))
        ctig = ctig[ctig %NN% dt2gr(rtmp[(seed)]) > 9]
        if(length(ctig)==0) return(NULL)
        mate = GenomicRanges::reduce(gr.sum.strand(gr.flipstrand(dt2gr(rtmp[!(seed)]))) %Q% (score > 0))
        mate = mate[mate %NN% gr.flipstrand(dt2gr(rtmp[!(seed)])) > 9]
        seed = mate %NN% ctig
        ctig = as.data.table(GenomicRanges::reduce(c(gr.fix(ctig, mate), gr.fix(mate[seed > 0], ctig))))
        ctig[, read.cov := dt2gr(ctig) %NN% dt2gr(rtmp[(seed)])] ##rtmp[, sum(seed)]]
        ctig = rbind(ctig, as.data.table(mate[seed==0])[, read.cov := rtmp[, sum(!seed)]], use.names=T, fill=T)
        ctig[, map60.cov := read.cov]
        ctig[, marker := c(0, cumsum(width))[1:.N]]
        ctig[, me := sum(width) - cumsum(width)]
        ctig[, cigar := paste0(ifelse(marker > 0, paste0(marker, "S"), ""), width, "M", ifelse(me > 0, paste0(me, "S"), ""))]
        ctig[strand=="-", cigar := paste0(ifelse(me > 0, paste0(me, "S"), ""), width, "M", ifelse(marker > 0, paste0(marker, "S"), ""))]
        ctig[, peak := gr.string(gr.stripstrand(win[,c()]))]
        ctig[, qname := i]
        return(ctig[, !c("marker", "me")])
    }), use.names=TRUE, fill=TRUE)

    if(nrow(out.dt)){
        out.dt[, track := t]
        out.dt[, leix := li$leix]
        out.dt[, qname := paste(li$sample, leix, track, qname, sep="_")]
        out.dt[, c_type := "human"]
        out.dt[, c_spec := "human"]
        out.dt$mapq = 60
        out.dt$seq = ""
    } else{
        out.dt$qname = character()
        out.dt$track = character()
        out.dt$leix = character()
        out.dt$c_type = character()
        out.dt$c_spec = character()
        out.dt$mapq = integer()
        out.dt$seq = character()
    }

    ##res = caller(li, out.dt[is.dup(qname)], ref_dir=ref_dir, ref_obj=ref_obj,
    res = caller(li, out.dt[is.dup(qname)], human = human,
                return.contigs = return.contigs)
    return(res)
}

#' transform
#'
#' this function allows conversion between data.table and gr.sum
#' (used for convenience with a data.table by= argument)
#' @param seq character seqnames
#' @param s integer start
#' @param e integer end
transform = function(seq, s, e){
    gr = GRanges(seq, IRanges(s, e))
    dt = gr2dt(gr.sum(gr))
    dt$seqnames = as.character(dt$seqnames)
    return(dt)
}


## #' .sample.spec
## #'
## #' loads reads and mates for a single sample (tumor or normal)
## #' @param le GRanges or data.table of loose ends
## #' @param bam path to BAM file
## #' @param pad integer width of padding to add around loose ends
## #' @param verbose optional, default=FALSE
## .sample.spec = function(le, bam, pad, verbose=FALSE){
##     if(verbose) message(paste("loading reads from", bam))
##     if(!inherits(le, "GRanges")) le = dt2gr(le)
##     sl = c(seqlengths(BamFile(bam)), setNames(1, "*"))
##     has.chr = any(grepl("chr", seqnames(seqinfo(BamFile(bam)))))
##     if(has.chr) { w = gr.chr(le) + pad
##     } else w = le + pad
##     reads = as.data.table(unlist(bamUtils::read.bam(bam, gUtils::gr.reduce(gr.stripstrand(w)), pairs.grl=T, isDuplicate=NA, isPaired=TRUE, tag="SA")))
##     splits = reads[!is.na(SA)]
##     if(nrow(splits) > 0){
##         splits$SA = as.character(splits$SA)
##         splwin = dunlist(strsplit(splits$SA, ";"))
##         spl = unlist(lapply(strsplit(splwin$V1, ","), function(w) paste(w[1], w[2], sep=":")))
##         spl = GRanges(spl)
##         spl$qname = splits[as.integer(splwin$listid)]$qname
##         splitsides = as.data.table(unlist(read.bam(bam, gUtils::gr.reduce(spl+150)-150, pairs.grl=T, isDuplicate=NA, tag="SA")) %Q% (qname %in% spl$qname))[order(mrnm, mpos)][!duplicated(paste(seqnames, start, qname, seq))]
##         reads = rbind(reads, splitsides, fill=T, use.names=TRUE)
##     }
##     reads[, unpmate := bamflag(flag)[, "hasUnmappedMate"]==1]
##     reads[, isunp := start == 1 & is.na(seq)]
##     reads[, unp := any(unpmate) & any(isunp), by=qname]
##     reads[(unp), ":="(start = ifelse(isunp, start[unpmate], start), end = ifelse(isunp, end[unpmate], end)), by=qname]
##     reads[, missing := any(is.na(seq)), by=qname]
##     mw = reads[!is.na(mrnm) & !is.na(seq)]
##     stopifnot(!is.na(mw$mrnm))
##     mw[, ":="(seqnames = mrnm, start = ifelse(is.na(mpos), start, mpos), end = ifelse(is.na(mpos), end, mpos))]
##     reads = reads[!is.na(seq)]
##     ##    fixmr = reads[qname %in% mw$qname, setNames(mrnm, qname)]
##     ##    mw[, seqnames := fixmr[qname]]
##     mw = dt2gr(mw[,!"strand"], seqlengths=sl)
##     mw[width(mw) == 0] = mw[width(mw) == 0] + 1
##     mate.wins = gUtils::gr.reduce(mw+150)-150
##     reads = reads[, !c("unpmate", "isunp", "unp", "SA")]
##     if(length(mate.wins) > 0){
##         if(verbose) message(paste("loading mates from", length(mate.wins), "windows"))
##         m.mate.wins = mate.wins
##         ## this is potentially slow for a large number of mate wins
##         mates = rbindlist(lapply(seq(1, length(m.mate.wins), 100), function(i){
##             mi = read.bam(bam, gUtils::gr.reduce(gr.stripstrand(m.mate.wins[i:min(length(m.mate.wins), i+99)])), pairs.grl=F, isDuplicate=NA)
##             mi = as.data.table(mi[mi$qname %in% reads$qname])
##             gc()
##             return(mi)
##         }), fill=TRUE, use.names=TRUE)
##         gc()
##         if(verbose) message("mates loaded")
##         rpair = rbind(reads, mates, fill=TRUE, use.names=TRUE)[!duplicated(paste(qname, flag)),]; rpair$MQ = NULL
##     } else {
##         rpair = reads[!duplicated(paste(qname, flag)),]; rpair$MQ = NULL
##     }
##     rpair[, R1 := bamflag(flag)[, "isFirstMateRead"]==1]
##     rpair[, R2 := bamflag(flag)[, "isSecondMateRead"]==1]
##     rpair[, paired := any(R1) & any(R2), by=qname]
##     if(verbose) message(ifelse(all(rpair$paired), "Found All Mates!!", "Some mates still missing - perhaps BAM was deduplicated"))
##     rpair[, MQ := rev(mapq), by=qname]
##     rpair[, count := .N, by = qname]
##     rpair[count == 0, MQ := 0]
##     ## rpair[, MQ := MQ * as.integer(.N>1), by=qname]
##     ##    flip = rpair$strand == "-"
##     flip = bamflag(rpair$flag)[, "isMinusStrand"] == 1 | rpair$strand == "-"
##     ##flip = bamflag(rpair$flag)[, "isMinusStrand"] == 1
##     rpair[flip, seq := as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seq)))]
##     ##    rpair = rpair[!duplicated(paste(qname, R1, seq))]
##     rpair = rpair[rev(order(nchar(seq)))][!duplicated(paste(qname, R1))]
##     reads = dt2gr(rpair)
##     gc()
##     return(reads)
## }

#' @name .sample.spec2
#' @title .sample.spec2
#'
#' @description
#' loads reads and mates for a single sample (tumor or normal)
#' assumes that BAM has already been filtered and avoids slow lapply
#' 
#' @param bam path to BAM file
#' @param chrsub (logical) substitute chr header? default TRUE
#' @param verbose optional, default=FALSE
.sample.spec2 = function(bam,
                         chrsub = TRUE,
                         verbose = FALSE) {
    if (verbose) {
        message(paste("loading reads from", bam))
    }

    ## load all sequences from BAM
    ## this assumes that this BAM has been pre-filtered to only include reads in relevant windows
    ## and their mates
    all.reads.grl = bamUtils::read.bam(bam, all = TRUE,
                                       pairs.grl = TRUE, ## return GRangesList with read pairs
                                       isDuplicate=NA, ## load all reads, regardless of if duplicated
                                       isPaired=TRUE,
                                       tag="SA") ## indicate split alignments
    reads = as.data.table(unlist(all.reads.grl))
    ## splits = reads[!is.na(SA)]
    ## if(nrow(splits) > 0){
    ##     splits$SA = as.character(splits$SA)
    ##     ## grab the windows into which the reads are split
    ##     splwin = dunlist(strsplit(splits$SA, ";"))
    ##     spl = unlist(lapply(strsplit(splwin$V1, ","), function(w) paste(w[1], w[2], sep=":")))
    ##     spl = GRanges(spl)
    ##     ## get the other side of the read with matching qname from the BAM file
    ##     spl$qname = splits[as.integer(splwin$listid)]$qname
    ##     splitsides = as.data.table(unlist(read.bam(bam, gUtils::gr.reduce(spl+150)-150, pairs.grl=T, isDuplicate=NA, tag="SA")) %Q% (qname %in% spl$qname))[order(mrnm, mpos)][!duplicated(paste(seqnames, start, qname, seq))]
    ##     reads = rbind(reads, splitsides, fill=T, use.names=TRUE)
    ## }
    reads[, unpmate := bamflag(flag)[, "hasUnmappedMate"]==1]
    reads[, isunp := start == 1 & is.na(seq)]
    reads[, unp := any(unpmate) & any(isunp), by=qname]
    ## for reads with that are unmapped
    ## set the start and end to start/end of its mate
    reads[(unp), ":="(start = ifelse(isunp, start[unpmate], start),
                      end = ifelse(isunp, end[unpmate], end)),
          by=qname]
    ## a missing read is an qname in which one of the sequences is NA
    reads[, missing := any(is.na(seq)), by=qname]
    reads = reads[!is.na(seq)]
    ## choose non-duplicated reads and designate one R1 and the other R2
    rpair = reads[!duplicated(paste(qname, flag)),];
    rpair$MQ = NULL
    rpair[, R1 := bamflag(flag)[, "isFirstMateRead"]==1]
    rpair[, R2 := bamflag(flag)[, "isSecondMateRead"]==1]
    rpair[, paired := any(R1) & any(R2), by=qname]
    reads = reads[, !c("unpmate", "isunp", "unp", "SA")]
    if(verbose) {
        message(ifelse(all(rpair$paired),
                       "Found All Mates!!",
                       "Some mates still missing - perhaps BAM was deduplicated"))
    }
    rpair[, MQ := rev(mapq), by=qname]
    rpair[, count := .N, by = qname]
    rpair[count == 0, MQ := 0]
    flip = bamflag(rpair$flag)[, "isMinusStrand"] == 1 | rpair$strand == "-"
    rpair[flip, seq := as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seq)))]
    rpair = rpair[rev(order(nchar(seq)))][!duplicated(paste(qname, R1))]
    reads = dt2gr(rpair)

    if (chrsub) {
        return(gr.nochr(reads))
    }
    gc()
    return(reads)
}
#' .realign
#'
#' @description
#' realigns reads and their mates single-end to attach individual MAPQs
#' 
#' @param reads (GRanges or data.table) reads from BAM file with $seq in original reading frame
#' @param ref (BWA object from RSeqLib) reference genome for realignment
#' @param gg (gGraph) sample used to identify sequences fitted in graph, default=NULL
#' @param filter (logical) return loose reads pairs only? default TRUE
#' @param chunksize (numeric) perform realignment in chunks to help with memory. default 1e6
#' @param verbose optional, default=FALSE
.realign = function(reads, ref, gg=NULL, filter=TRUE, verbose=F, chunksize = 5e5){
    if(!is.null(gg)){
        seqs = unique(seqnames(gg$nodes[!is.na(cn)]$gr))
    } else {
        seqs = c(1:22, "X", "Y")
        if(any(grepl("chr", seqlevels(ref)))) seqs = gr.chr(seqs)
    }

    ## realign ALL reads if not filtering for loose pairs
    ## otherwise, realign just the low MAPQ read in the pair
    if(filter){
        qni = as.data.table(reads)[is.na(mapq) | mapq<50 | is.na(MQ), (unique(qname))]
    } else{
        qni = unique(reads$qname)
    }

    ## start realignment (in chunks)
    redo = reads$qname %in% qni
    if(verbose) {
        message(paste("realigning", sum(redo), "reads"))
    }
    
    ## do realignment in chunks to help with memory
    seqs.for.realignment = setNames(reads$seq, 1:length(reads))
    ix.for.realignment = which(redo)
    nchunks = ceiling(sum(redo) / chunksize)
    realn = data.table()
    for (chunk in 1:nchunks) {
        starti = (chunk - 1) * chunksize + 1
        endi = pmin(sum(redo), (chunk * chunksize))
        if (verbose) {
            message("Chunk start: ", starti)
            message("Chunk end: ", endi)
            message("Chunk ", chunk, " of ", nchunks)
        }
        realni = ref[seqs.for.realignment[ix.for.realignment[starti:endi]]]
        gc()
        realn = rbind(realn, as.data.table(realni), use.names = TRUE, fill = TRUE)
        ## remove unneeded variables to help with memory footprint
        rm(realni)
        gc()
    }
    
    ## recast as GRanges since legacy code expects that
    realn = dt2gr(realn)
    ## realn = ref[setNames(reads$seq, 1:length(reads))[redo]]
    realn$mapq = as.integer(realn$mapq)
    realn$query.id = as.integer(realn$qname); realn$qname = reads[realn$query.id]$qname
    values(realn) = cbind(values(realn), values(reads[realn$query.id, !(colnames(values(reads)) %in% colnames(values(realn)))]))

    ## keep track of the unaligned reads that didn't appear in the GRanges from BWA
    uix = ! ( (1:length(reads))[redo] %in% realn$query.id )
    uix = (1:length(reads))[redo][uix]
    unaln = as.data.table(reads[uix])
    unaln[, ":="(
        seqnames = "*",
        start = 1,
        end = 0,
        flag = ifelse(R1, 69, 113),
        mapq = 0,
        query.id = uix)]

    ## final reads are the realigned reads in new genomic coordanates + unaligned reads
    ## in addition, the original sequence + realignment mapqs are stored
    realn = rbind(as.data.table(realn), unaln, fill=T, use.names=TRUE)
    realn$reading.frame = reads[realn$query.id]$seq
    realn$mapq = as.integer(realn$mapq)
    realn[, mapq := ifelse(!(seqnames %in% seqs), as.integer(0), mapq), by=query.id]
    realn = realn[rev(order(mapq))][!duplicated(query.id), ]
    gc()
    realn[, MQ := ifelse(rep(.N, .N)==1, as.integer(NA), c(mapq[-1], mapq[1])), by=qname]
    cols = c(colnames(realn)[colnames(realn) %in% colnames(as.data.table(reads))], "reading.frame", "AS")
    realn = realn[, cols, with=F]

    ## annotate whether the read belongs to a loose read pair
    lqn = realn[mapq>50 & (is.na(MQ) | MQ < 1), qname]
    if(filter){
        realn = realn[qname %in% lqn,]
        realn[, loose.pair := TRUE]
    } else realn[, loose.pair := qname %in% lqn]
    realn[, high.mate := mapq>50 & (is.na(MQ) | MQ < 1)]
    gc()
    return(realn)
}


## #' loose.reads
## #'
## #' loads and reads and their mates from given windows and realigns single end for read-specific MAPQ
## #' @param le GRanges or data.table windows around which to sesarch for reads
## #' @param tbam character path to tumor BAM file
## #' @param pad optional, integer padding around le windows, default=25e3
## #' @param nbam optional, character path to normal BAM file, default=NULL
## #' @param ref optional, BWA object of reference genome for realignment or path to fasta, default=package/extdata/hg19_looseends/human_g1k_v37_decoy.fasta
## #' @param filter optional, logical filter=TRUE returns loose read pairs only, filter=FALSE returns all read pairs annotated with logical $loose column, default=T
## #' @param gg optional, gGraph corresponding to sample, used to identify sequences fitted in graph, default=NULL
## #' @param verbose optional, default=FALSE
## #' @export
## loose.reads = function(le, tbam, pad=25e3, nbam=NA, ref=system.file('extdata', 'hg19_looseends', 'human_g1k_v37_decoy.fasta', package='loosends'), filter=TRUE, gg=NULL, verbose=FALSE){
##     le = copy(le)
##     if(inherits(ref, "character")){
##         if(!file.exists(ref)) stop("Provide reference BWA object or path to reference fasta for loose.reads")
##         ref = RSeqLib::BWA(fasta=ref)
##     }
##     if(is.na(nbam)) nbam = NULL
##     id = le[1]$sample
##     treads = .sample.spec(copy(le), tbam, pad, verbose=verbose)
##     realn = .realign(treads, ref, filter=filter, gg=gg, verbose=verbose)
##     realn$sample = id

##     if(!is.null(nbam)){
##         nreads = .sample.spec(copy(le), nbam, pad, verbose=verbose)
##         nrealn = .realign(nreads, ref, filter=filter, gg=gg, verbose=verbose)
##         nrealn$sample = paste0(id, "N")
##         gc()
##         return(rbind(realn, nrealn, fill=TRUE, use.names=TRUE))
##     }
##     gc()
##     return(realn)
## }


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
    

#' @name grab_ref_obj
#' @title grab_ref_obj
#'
#' @param ref.dir (character) path to reference directory
#'
#' @return list with names rep, human, polyA, microbe
grab_ref_obj = function(ref.dir = NA_character_) {

    if (!dir.exists(ref.dir)) {
        stop("Supplied directory does not exist")
    }

    human = BWA(fasta=paste(ref.dir, "human_g1k_v37_decoy.fasta", sep="/"))
    rep = BWA(fasta=paste(ref.dir, "mskilab_combined_TraFicv8-3_satellites.fa", sep="/"))
    polyA = BWA(fasta=paste(ref.dir, "PolyA.fa", sep="/"))
    microbe = BWA(fasta=paste(ref.dir, "human_g1k_v37.withviral.fasta", sep="/"),
                  keep_sec_with_frac_of_primary_score=0.2)

    ref.obj = list(human=human, rep=rep, polyA=polyA, microbe=microbe)
    return(ref.obj)
}

#' @name call_loose_end
#' @title call_loose_end
#'
#' @description
#'
#' This function calls a loose end given reads from near to the loose end
#'
#' @param li (data.table) data.table coercible to GRanges, must have metatdata $leix and $sample
#' @param ri (data.table) loose reads table prepped (e.g. from prep_loose_reads)
#' @param concat.bwa (BWA object from RSeqLib)
#' @param human.bwa (BWA object from RSeqLib)
#' @param pad (numeric) window size for local assembly, default 1 kbp
#' @param mix.tn (logical) mix tumor/normal reads before assembly
#' @param max.reads (logical) maximum number of reads before downsampling (default 5000)
#' @param verbose
#'
#' @return list with entries
#' - all.contigs (data.table)
#' - wide.contigs (data.table or NULL)
#' - call (data.table)
call_loose_end = function(li, ri,
                          concat.bwa = NULL, ## BWA object with concatenated reference
                          human.bwa = NULL, ## BWA object with just human reference
                          concat.fn = "~/git/loosends/inst/extdata/hg19_loosends/concatenated_references_deduped.fasta",
                          pad = 5e3,
                          mix.tn = TRUE,
                          max.reads = 2500,
                          outdir = "./",
                          minimap = FALSE,
                          verbose = FALSE) {
    
    ## human = ref_obj$human
    ## rep = ref_obj$rep
    ## polyA = ref_obj$polyA
    ## microbe = ref_obj$microbe
    
    uannot = readRDS(system.file('extdata', '101.unmappable.annotations.rds', package='loosends'))

    id = li$sample

    if (!dir.exists(outdir)) {
        if (verbose) message("Creating output directory")
        dir.create(outdir, recursive = TRUE)
    }

    if (verbose) {
        message("Generating call")
    }

    .build.tigs = function(ri, pp, id, leix, verbose = TRUE) {
        if(verbose) message("assembling contigs")
        out.dt = rbindlist(lapply(1:length(pp), function(i){
            win = pp[i]
            build.from.win(win, ri, verbose = verbose)
        }), use.names=TRUE, fill=TRUE)
        if("seed" %in% colnames(ri)) ri$seed = NULL
        out.dt$sample = rep(id, nrow(out.dt))
        out.dt$leix = rep(leix, nrow(out.dt))
        out.dt$somatic = rep(somatic, nrow(out.dt))
        out.dt = out.dt[!is.na(seq)]

        if(verbose) message(paste("assembled", nrow(out.dt), "contigs across all tracks"))
        if(nrow(out.dt)){
            out.dt[, qname := 1:.N, by=.(track, sample, leix)]
            out.dt[, qname := paste(sample, leix, track, qname, sep="_")]
        } else{
            out.dt$qname = character()
        }

        ## save contigs as FASTA
        if (out.dt[, .N]) {
            unique.tigs = unique(out.dt, by = "qname")
            tigs.vec = setNames(object = unique.tigs[, seq], nm = unique.tigs[, qname])
            tigs.strings = Biostrings::DNAStringSet(tigs.vec)
            Biostrings::writeXStringSet(tigs.strings, paste0(outdir, "/", "contigs.fasta"))
        }

        if (verbose) {
            message("Aligning contigs to concatenated reference")
        }

        if (minimap) {
            if (verbose) message("Using minimap for contig alignment")

            contigs.aln.bam.fn = paste0(outdir, "/", "contigs.aln.bam")
            contigs.aln.sorted.bam.fn = paste0(outdir, "/", "contigs.aln.sorted.bam")
            contigs.aln.sam.fn = paste0(outdir, "/", "contigs.aln.sam")

            if (verbose) message("starting short read single end alignment")
            cmd = paste("minimap2 -ax sr --secondary=yes",
                        concat.fn,
                        paste0(outdir, "/", "contigs.fasta"),
                        ">",
                        contigs.aln.sam.fn)

            sys.res = system(cmd)
            if (sys.res) {
                file.remove(contigs.aln.sam.fn)
                stop("Error!")
            }

            cmd = paste("samtools view -Sb", contigs.aln.sam.fn, ">", contigs.aln.bam.fn)
            sys.res = system(cmd)
            
            if (sys.res) {
                file.remove(contigs.aln.bam.fn)
                stop("Error!")
            }

            cmd = paste("samtools sort", contigs.aln.bam.fn, ">", contigs.aln.sorted.bam.fn, ";",
                        "samtools index", contigs.aln.sorted.bam.fn)
            sys.res = system(cmd)
            
            if (sys.res) {
                file.remove(contigs.aln.bam.sorted.fn)
                stop("Error!")
            }

            ## read BAM from single end alignment
            aln.gr = read.bam(bam = contigs.aln.sorted.bam.fn,
                              all = TRUE,
                              pairs.grl = FALSE,
                              tag="AS",
                              isPaired = NA)
            ralns.og = aln.gr
            ## don't include unaligned reads?
            if (length(ralns.og)) {
                ralns.og = ralns.og %Q% (!is.na(qname))
            }

            ## garbage collect afterwards to avoid memory issues??
            ## hack hack
            gc()
        } else {
            ralns.og = concat.bwa[out.dt$seq]
            if (length(ralns.og)) {
                values(ralns.og)[, "qname"] = out.dt[as.integer(ralns.og$qname), qname]
            }
        }

        ## TODO: store this as a package global variable
        keepseq = seqnames(seqinfo(concat.bwa))[c(1:87, 148:6398)]
        rep.seqs = seqnames(seqinfo(concat.bwa))[1:60]
        poly.seqs = seqnames(seqinfo(concat.bwa))[61]
        human.seqs = seqnames(seqinfo(concat.bwa))[62:86]
        viral.seqs = seqnames(seqinfo(concat.bwa))[148:6398]

        
        keep.cols = c("cigar", "flag", "mapq", "AS", "qname")

        if (length(ralns.og)) {
            ralns.dt = as.data.table(ralns.og[, keep.cols]  %Q% (as.character(seqnames(ralns.og)) %in% keepseq))

            ## transfer original values based on qname
            setnames(ralns.dt, old = "qname", new = "qname.og")
            out.dt.cols = setdiff(colnames(out.dt), colnames(ralns.dt))
            ralns.dt = cbind(ralns.dt, out.dt[match(ralns.dt$qname.og, qname), ..out.dt.cols])

            ralns.dt[as.character(seqnames) %in% rep.seqs, c_type := "rep"]
            ralns.dt[as.character(seqnames) %in% poly.seqs, c_type := "polyA"]
            ralns.dt[as.character(seqnames) %in% human.seqs, c_type := "human"]
            ralns.dt[as.character(seqnames) %in% viral.seqs, c_type := "viral"]
            ralns.dt[, c_spec := c_type]
            ralns.dt[c_type == "rep", c_spec := dunlist(strsplit(as.character(seqnames), "#", fixed=T))[rev(!duplicated(rev(listid)))][order(as.integer(listid)), V1]]
            calns = ralns.dt[!(is.na(cigar) | is.na(AS) | is.na(flag) | is.na(mapq))]
        } else {
            calns = as.data.table(ralns.og)
        }

        ## if(verbose) message("aligning contigs to 4 references")
        ## ureps = rep[out.dt$seq]
        ## preps = polyA[out.dt$seq]
        ## mreps = microbe[out.dt$seq]
        ## hreps = human[out.dt$seq]
        ## virals = seqlevels(microbe)[85:length(seqlevels(microbe))]
        
        ## if(!(length(ureps))){
        ##     ureps$cigar = character()
        ##     ureps$mapq = character()
        ##     ureps$AS = integer()
        ##     ureps$flag = character()
        ## }
        ## if(!(length(preps))){
        ##     preps$cigar = character()
        ##     preps$mapq = character()
        ##     preps$AS = integer()
        ##     preps$flag = character()
        ## }
        ## if(!(length(hreps))){
        ##     hreps$cigar = character()
        ##     hreps$mapq = character()
        ##     hreps$AS = integer()
        ##     hreps$flag = character()
        ## }
        ## if(!(length(mreps))){
        ##     mreps$cigar = character()
        ##     mreps$mapq = character()
        ##     mreps$AS = integer()
        ##     mreps$flag = character()
        ## }
        ## keep.cols = c("cigar", "flag", "mapq", "AS")
        ## values = rbind(out.dt[as.integer(ureps$qname)], out.dt[as.integer(preps$qname)], out.dt[as.integer(hreps$qname)], fill=T, use.names=T)
        ## ralns = rbind(as.data.table(ureps[, keep.cols]), as.data.table(preps[, keep.cols]), as.data.table(hreps[, keep.cols]))
        ## ralns = cbind(ralns, values)
        ## if(nrow(out.dt)){
        ##     rch = cgChain(ralns)
        ##     good.ids = c(as.character(seqnames(gaps(rch$x) %Q% (strand=="+"))), setdiff(out.dt$qname, ralns$qname))
        ##     mreps$query.id  = as.integer(mreps$qname)
        ##     mreps$qname = out.dt[mreps$query.id, qname]
        ##     vreps = mreps %Q% (seqnames %in% virals) %Q% (qname %in% good.ids)
        ##     if(verbose & length(vreps)) message("adding viral alignments")
        ##     valns = cbind(as.data.table(vreps[, keep.cols]), out.dt[vreps$query.id])
        ##     calns = rbind(ralns, valns)
        ##     calns[, c_type := c(rep("rep", length(ureps)), rep("polyA", length(preps)), rep("human", length(hreps)), rep("viral", length(vreps)))]
        ##     calns[, c_spec := c_type]
        ##     calns[c_type == "rep", c_spec := dunlist(strsplit(as.character(seqnames), "#", fixed=T))[rev(!duplicated(rev(listid)))][order(as.integer(listid)), V1]]
        ## } else{
        ##     calns = ralns
        ## }
        calns$somatic = rep("somatic", nrow(calns))
        calns$mapq = as.integer(calns$mapq)

        ## fill in dummy values
        calns[, ":="(read.cov = NA, map60.cov = NA, tumor.frac = NA,
                     tumor.support = NA, normal.support = NA, total.support = NA)]


        if(verbose) message("parsing contigs for telomeric matches")
        all.contigs = copy(calns)
        if (nrow(all.contigs)) {
            all.contigs[c_type == "human", c_spec := ifelse(seqnames %in% c(1:22, "X", "Y"), "human", "unassembled")]
            all.contigs$cmer = munch(all.contigs, eighteenmer('c'))
            all.contigs$gmer = munch(all.contigs, eighteenmer('g'))

            if(verbose) message("parsing alignments for unmappable repeat overlaps")
            coord.calls = copy(all.contigs[c_type=="human"])[, !c("c_type", "c_spec")]
            ov = gr.findoverlaps(dt2gr(coord.calls), uannot)
            if(length(ov)){
                subj = ov$subject.id
                strand(ov) = coord.calls[ov$query.id, strand]
                values(ov) = values(dt2gr(coord.calls)[ov$query.id])
                ov$c_type = "rep"
                ov$c_spec = uannot[subj]$c_spec
                coord.calls = as.data.table(ov)
                return(rbind(all.contigs, coord.calls, fill=T, use.names=T))
            }
        }
        return(all.contigs)
    }

    if (verbose) {
        message("Building contigs")
    }

    somatic = as.logical(nrow(ri[grepl("control", track)]))
    wholseed = dt2gr(li)[,c()] + pad ##dt2gr(li)[,c()]+1e3
    if (verbose) {
        message("Using window size: ", pad, " bp")
    }
    pp = gr.stripstrand(gr.tile(wholseed, 200) %Q% (width == 200))

    ## check if contigs should be assembled mixing tumor and normal reads
    ri.input = copy(ri)
    if (mix.tn) {
        if (verbose) {
            message("Mixing tumor/normal reads before assembly")
        }
        ri.input$track = NULL ## set track to null so that assembly is sample-agnostic
    }
    ## check if reads need to be downsampled before assembly
    if (ri.input[, .N]) {
        n.qn = length(unique(ri.input[, qname]))
        if (verbose) {
            message("Number of unique read qnames for assembly: ", n.qn)
        }
        if (n.qn > max.reads) {
            if (verbose) {
                message("Downsampling reads to max # qnames: ", max.reads)
            }
            keep.qn = sample(unique(ri.input[, qname]), max.reads, replace = FALSE)
            ri.input = ri.input[qname %in% keep.qn,]
        }
    }
    all.contigs = .build.tigs(ri.input, pp, li$sample, li$leix, verbose = verbose)

    ## also save the filtered contigs
    if (verbose) {
        message("Generating call")
    }

    if (mix.tn) {
        res = caller(li, all.contigs, human = human.bwa, 
                     return.contigs = TRUE, pad = 5e3)
    } else {
        ## if not mixed, call should only be generated from tumor contigs
        res = caller(li, all.contigs[grepl('sample', track),],
                     human = human.bwa,
                     return.contigs = TRUE)
    }


    ## add junction annotation to cal
    call = res$call
    filtered.contigs = res$contigs
    ## HACK HACK
    bowtie = max(ri[, mapq], na.rm = TRUE) < 60 ## you know it's bowtie if mapqs only go up to 44
    if (bowtie) {
        recall.res = read.based(li, ri, pad = pad, human = human.bwa, return.contigs = TRUE, mapq.thresh = 30)
    } else {
        recall.res = read.based(li, ri, pad = pad, human = human.bwa, return.contigs = TRUE)
    }
    disc = recall.res$call
    recall = (!call$missedj & disc$missedj) | (!call$complex & disc$complex)
    call[, missedj := missedj | disc$missedj]
    call[, complex := complex | disc$complex]
    call[junction == "" & disc$junction != "", junction := disc$junction]
    call[junction != "" & disc$junction != "", junction := {
        cgr = parse.gr(junction)
        dgr = parse.gr(disc$junction)
        s = gr.reduce(cgr[1], dgr[1], ignore.strand=F)
        m = gr.reduce(cgr[-1], dgr[-1], ignore.strand=F)
        m = gr.start(GenomicRanges::reduce(m + 200) - 200, ignore.strand=F)
        paste(c(gr.string(s), gr.string(m)), collapse=" | ")
    }]
    call[, mystery := !missedj & !complex & mate.mappable & seed.mappable & !insertion]
    if(any(recall)){
        call[recall]$call = update.call(call[recall])
        call[(missedj), seed.mappable := TRUE]
        call[(missedj), mate.mappable := TRUE]
        filtered.contigs = grbind(filtered.contigs, recall.res$contigs, fill = TRUE, use.names = TRUE)
    }
    ## if(call$mystery){
    ##     if(verbose) message("mystery: repeating assembly at larger intervals")
    ##     pp = gr.stripstrand((gr.tile(wholseed-250, 500)+250) %Q% (width == 1e3))
    ##     ## changed from ri to ri.input for consistency
    ##     wide.contigs = .build.tigs(ri.input, pp, li$sample, li$leix, verbose = TRUE)
    ##     if (mix.tn) {
    ##         call = caller(li, wide.contigs, ref_obj=ref_obj, pad = 5e3)
    ##     } else {
    ##         call = caller(li, wide.contigs[grepl('sample', track),], ref_obj = ref_obj, pad = 5e3)
    ##     }
    ## }

    call[, category := ifelse(complex, "complex rearrangement",
                       ifelse(missedj, "missed junction",
                       ifelse(mystery | grepl("mystery", mate.repeats) | (seed.mappable & mate.mappable),
                              "mystery",
                              paste("type", as.integer(!seed.mappable) + as.integer(!mate.mappable),
                                    "loose end"))))]

    ## update call to set category used in paper (type 0 loose ends)
    call[, paper_category := ifelse(category == "complex rearrangement" | category == "missed junction",
                                    "type 0 loose end",
                                    category)]
    
    call[category=="mystery", mystery := TRUE]
    
    call[,  ":="(
        leix = li[, leix],
        loose.end = li[, paste0(seqnames, ":", start, strand)],
        sample = id)]

    if (verbose) {
        message("category: ", call$category)
    }

    res = list(all.contigs = all.contigs,
               call = call,
               filtered.contigs = as.data.table(filtered.contigs))
    ## if (exists('wide.contigs')) {
    ##     if (inherits(wide.contigs, "data.table")) {
    ##         res$wide.contigs = wide.contigs
    ##     }
    ## }
    return(res)
}
