## #' contig.support
## #'
## #' quantifies preferential read support for contig over reference
## #' @param reads GRanges of seed reads and mates
## #' @param contig data.table of reference alignments for a single contig
## #' @param ref optional, character of reference sequence surrounding loose end, default=NULL
## #' @param min.bases optional, minimum aligned bases from a read to count for support, default=20
## #' @param min.aligned.frac optional, minimum fraction of aligned bases from a read to count for support, default=0.95
## contig.support.legacy = function(reads, contig, ref = NULL, min.bases = 20, min.aligned.frac = 0.95) {
##     cg.contig = gChain::cgChain(contig)
##     if (length(reads) == 0) 
##         stop("reads must be non empty GRanges with $qname, $cigar, $seq, and $flag fields")
##     if (length(contig) == 0) 
##         stop("contig must be non empty GRanges with $qname, $cigar and $seq fields")
##     seq = unique(gr2dt(contig), by = c("qname"))[, structure(as.character(seq), 
##                                                              names = as.character(qname))]
##     bwa.contig = RSeqLib::BWA(seq = seq)
##     chunks = gChain::links(cg.contig)$x
##     strand(chunks) = "+"
##     chunks = disjoin(chunks)
##     if(!("R1" %in% colnames(values(reads))))
##         reads$R1 = bamflag(reads$flag)[, "isFirstMateRead"] > 0
##     if (is.null(reads$AS)) {
##         warning("AS not provided in reads, may want to consider using tag = \"AS\" argument to read.bam or provide a ref sequence to provide additional specificity to the contig support")
##         reads$AS = 0
##     }
##     nix = as.logical(strand(reads) == "-")
##     reads$seq[nix] = reverseComplement(DNAStringSet(reads$seq[nix]))
##     reads[!reads$R1] = gr.flipstrand(reads[!reads$R1])
##     reads$seq[!reads$R1] = reverseComplement(DNAStringSet(reads$seq[!reads$R1]))
##     reads = reads %Q% (!duplicated(paste(qname, R1)))
##     if (!is.null(ref)) {
##         bwa.ref = RSeqLib::BWA(seq = ref)
##         tmp = bwa.ref[reads$seq] %>% gr2dt
##         tmp$ix = as.numeric(as.character(tmp$qname))
##         tmp$R1 = reads$R1[tmp$ix]
##         tmp$qname = reads$qname[tmp$ix]
##         tmp = unique(tmp, by = c("qname", "R1"))
##         setkeyv(tmp, c("qname", "R1"))
##         if (nrow(tmp)) 
##             reads$AS = tmp[.(reads$qname, reads$R1), AS]
##     }
##     readsc = bwa.contig[reads$seq] %>% gr2dt
##     readsc$ix = as.integer(as.character(readsc$qname))
##     readsc$R1 = reads$R1[as.numeric(readsc$ix)]
##     readsc[, `:=`(nsplit, .N), by = .(qname, R1)]
##     readsc[, `:=`(aligned, countCigar(cigar)[, "M"])]
##     readsc[, `:=`(aligned.frac, aligned/qwidth[1]), by = .(qname, 
##                                                            R1)]
##     readsc$AS.og = reads$AS[readsc$ix]
##     readsc$AS.og[is.na(readsc$AS.og)] = 0
##     readsc$qname = reads$qname[readsc$ix]
##     ov = dt2gr(readsc) %*% chunks
##     strand(ov) = readsc$strand[ov$query.id]
##     ov$subject.id = paste0("chunk", ov$subject.id)
##     ovagg = dcast.data.table(ov %>% gr2dt, qname ~ subject.id, 
##                              value.var = "width", fun.aggregate = sum)
##     ovagg$nchunks = rowSums(ovagg[, -1] > min.bases)
##     rstats = gr2dt(ov)[, .(contig.id = unique(seqnames)[1], pos = sum(width[strand == 
##                                                                             "+"]), neg = sum(width[strand == "-"]), aligned.frac = min(aligned.frac), 
##                            num.contigs = length(unique(seqnames)), isize.contig = diff(range(c(start, 
##                                                                                                end))), qsplit = any(nsplit > 1), worse = any(AS.og > 
##                                                                                                                                              AS)), by = qname] %>% merge(ovagg, by = "qname")
##     keepq = rstats[nchunks > 1 & (pos == 0 | neg == 0) & aligned.frac > 
##                    min.aligned.frac & !worse & !qsplit & num.contigs == 
##                    1, ]
##     if (nrow(keepq) == 0) 
##         return(reads[c()])
##     keepq$aligned.frac = NULL
##     readsc = merge(readsc, keepq, by = "qname") %>% dt2gr
##     out = gChain::lift(cg.contig, readsc)
##     out[!out$R1] = gr.flipstrand(out[!out$R1])
##     out$col = ifelse(out$R1, "blue", "gray")
##     out
## }
