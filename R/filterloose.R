#' @importFrom stats pnorm t.test wilcox.test setNames

#' @name rel2abs
#' @title rel2abs
#'
#' rescales CN values from relative to "absolute" (i.e. per cancer cell copy) scale given purity and ploidy (forked from skitoosl)
#'
#' takes in gr with signal field "field"
#'
#' @param gr GRanges input with meta data field corresponding to mean relative copy "mean" in that interval
#' @param purity purity of sample
#' @param ploidy ploidy of sample
#' @param gamma gamma fit of solution (over-rides purity and ploidy)
#' @param beta beta fit of solution (over-rides purity and ploidy)
#' @param field meta data field in "gr" variable from which to extract signal, default "ratio"
#' @param field.ncn meta data field in "gr" variable from which to extract germline integer copy number, default "ncn", if doesn't exist, germline copy number is assumed to be zero
#' @param data_mean optionally provide a mean value to use in the transformation. Usually a mean value is computed from the input data, but in unique cases, where the input data does not represent the full set of data, then this value could be provided. For example, this is usefull when transforming SNV read counts, you can provide the ALT read count as the input that you want transformed, and provide the average count of ALT + REF as the data_mean
#' @param ncn.gr GRanges with the copy number values for normal samples (if the field.ncn is found in the input gr then the ncn.gr parameter is ignored). Notice the the input ncn.gr must contain the field specified by field.ncn
#' @param allele (logical) whether to return allelic CNs. If TRUE, assumes that the GRanges is "melted" and there are two identical ranges per SNP. Default FALSE.
#' @return
#' numeric vector of integer copy numbers
#' @export
rel2abs = function(gr, purity = NA, ploidy = NA, gamma = NA, beta = NA, field = 'ratio', field.ncn = 'ncn', data_mean = NA, ncn.gr = NA, allele = FALSE)
{
  mu = values(gr)[, field]
  mu[is.infinite(mu)] = NA
  w = as.numeric(width(gr))
  w[is.na(mu)] = NA
  sw = sum(w, na.rm = T)
  if (is.na(data_mean)){
      data_mean = sum(mu * w, na.rm = T) / sw
  }

  ncn = NA
  if (!is.null(field.ncn))
    if (field.ncn %in% names(values(gr)))
      ncn = values(gr)[, field.ncn]

  if (is.na(ncn)){
      if (!is.na(ncn.gr)){
          if (!inherits(ncn.gr, 'GRanges')){
              stop('ncn.gr must be of class GRanges, but ', class(GRanges), ' was provided.')
          }
          ncn = values(gr %$% ncn.gr[, field.ncn])[, field.ncn]
      } else {
      ncn = rep(2, length(mu))
      }
  }


  ploidy_normal = sum(w * ncn, na.rm = T) / sw  ## this will be = 2 if ncn is trivially 2

  if (allele) {
      y.bar = ploidy_normal * data_mean
      denom = purity * ploidy + ploidy_normal * (1 - purity)
      if (is.na(beta)) {
          beta = y.bar * purity / denom
      }
      if (is.na(gamma)) {
          gamma = (y.bar * (1 - purity)) / denom
      }
      return ((mu - gamma) / beta)
  }


  if (is.na(gamma))
    gamma = 2*(1-purity)/purity

  if (is.na(beta))
    beta = ((1-purity)*ploidy_normal + purity*ploidy) / (purity * data_mean)

  return(beta * mu - ncn * gamma / 2)
}

#' @name filter_loose
#' @title filter_loose
#'
#' @param jabba_rds (character) path to JaBbA
#' @param coverage (character) path to coverage file with gc-adjusted counts
#' @param norm.field (character) normal count field (reads.corrected)
#' @param mask (character) path to coverage mask
#' @param mask_pad (numeric) pad on distance to mask (in practice just for annotation purposes)
#' @param norm_thresh (numeric) 0.6 (from filter.loose)
#' @param min_size (numeric) minimum size in BP
#' @param subset (logical) subset loose ends to include only valid ones before returning
#' @param verbose (logical)
filter_loose = function(jabba_rds,
                      coverage = "/dev/null",
                      norm.field = 'reads.corrected',
                      mask = "~/projects/gGnome/files/zc_stash/tile.mask.rds",
                      mask_pad = 0,
                      norm_thresh = 0.6,
                      min_size = 0,
                      subset = FALSE,
                      verbose = TRUE) {

    if (!check_file(jabba_rds)) {
        stop("JaBbA file path not valid")
    }

    if (!check_file(coverage)) {
        if (verbose) {
            message("Coverage file not supplied!")
        }
        use.coverage = FALSE
    } else {
        if (verbose) {
            message("Reading normal coverage")
        }
        cov.gr = readRDS(coverage)
        cov.gr = gr.nochr(cov.gr)
        use.coverage = TRUE
    }

    if (use.coverage) {
        if (!(norm.field %in% names(values(cov.gr)))) {
            stop("supplied field not in coverage metadata: ", norm.field)
        }
        cov.gr = cov.gr %Q% (!is.infinite(values(cov.gr)[[norm.field]]))
    }

    ## pull loose ends from jabba_rds
    if (verbose) {
        message("Grabbing loose ends from JaBbA input")
    }
    le.gr = gg_loose_end(jabba_rds, return.type = 'GRanges')
    le.dt = as.data.table(le.gr)

    if (!nrow(le.dt)) {
        if (verbose) {
            message("No loose ends!")
        }
        return(le.dt)
    }
    
    le.gr$loose.index = 1:length(le.gr)

    ## check if loose end overlaps mask
    if (check_file(mask)) {

        if (verbose) {
            message("Checking overlaps with mask")
        }
        mask.gr = readRDS(mask)
        overlap.mask = le.gr %^% (gr.stripstrand(mask.gr) + mask_pad)
        le.dt[, mask := overlap.mask]

        ## also remove masked coverage bins!
        if (use.coverage) {
            cov.gr = cov.gr %Q% (!(cov.gr %^% gr.stripstrand(mask.gr)))
        }

        ## remove masked background bins
        if (use.bg) {
            bg.gr = bg.gr %Q% (!(bg.gr %^% gr.stripstrand(mask.gr)))
        }
        
    } else {
        le.dt[, mask := FALSE]
    }

    ## get node ids of fused and unfused loose ends
    if (verbose) {
        message("Identifying fused and unfused sides of each loose end")
    }

    ## grab gGraph
    gg = gG(jabba = jabba_rds)

    ## re2labs transformation
    if (use.coverage) {
        if (verbose) { message("rel2abs transforming normal coverage") }
        cov.gr$abscn = rel2abs(cov.gr, field = norm.field, purity = 1, ploidy = 2)
    }

    fused.unfused.dt = fused_unfused(le.dt, jabba_rds)

    ## get flanking normal coverage (borrowed from filter.loose)
    fused.melted.dt = melt.data.table(fused.unfused.dt, id.vars = "loose.index",
                                   measure.vars = c("fused.node.id", "unfused.node.id"),
                                   variable.name = "fused",
                                   value.name = "node.id")[!is.na(node.id),]
    
    segs.gr = gg$nodes$gr[fused.melted.dt$node.id, c()]
    names(segs.gr) = NULL
    segs.gr$loose.index = fused.melted.dt$loose.index
    segs.gr$fused = fused.melted.dt$fused == "fused.node.id"

    if (use.coverage) {

        ## use breakpoint coverage function?
        ## only keep coverage within 1e5 of loose end?
        sub.cov.gr = cov.gr %Q% (cov.gr %^% (le.gr + 1e5))
        ## sides.cov.dt = gr.findoverlaps(cov.gr, segs.gr, return.type = 'data.table')
        sides.cov.dt = gr.findoverlaps(sub.cov.gr, segs.gr, return.type = 'data.table')
        if (nrow(sides.cov.dt)) {

            ## use rel2abs transformation
            sides.cov.dt[, normal.cov := values(sub.cov.gr)[query.id, "abscn"]]
            ## use foreground (this should be always positive and can be logged)
            sides.cov.dt[, normal.fg := values(sub.cov.gr)[query.id, norm.field]]
            sides.cov.dt[, fused := values(segs.gr)$fused[subject.id]]
            sides.cov.dt[, loose.index := values(segs.gr)$loose.index[subject.id]]
        }

    }

    if (use.coverage) {
        ## compute mean coverage for fused and unfused sides

        if (sides.cov.dt[, .N]) {
            mean.cov.dt = sides.cov.dt[,
                                       .(mean.normal.cov = mean(.SD$normal.cov, na.rm = T)),
                                       by = .(loose.index, fused)]
            mean.cov.dt[, fused := ifelse(fused, 'fused', 'unfused')]
            mean.cov.dt = dcast.data.table(mean.cov.dt, loose.index ~ fused, value.var = "mean.normal.cov")

            ## also do KS-test
            sides.cov.dt[, both.sides := sum(.SD$fused == FALSE) > 3 & sum(.SD$fused == TRUE) > 3, by = loose.index]
            sides.cov.dt[, all.same := (length(unique(.SD$normal.fg[.SD$fused == TRUE])) < 2) |
                               (length(unique(.SD$normal.fg[.SD$fused == FALSE])) < 2)]
            if (sides.cov.dt[(both.sides), .N]) {
                ks.cov.dt = sides.cov.dt[(both.sides) & (!all.same),
                                         .(ks = wilcox.test(.SD$normal.fg[.SD$fused == TRUE],
                                                            .SD$normal.fg[.SD$fused == FALSE])$p.value,
                                           t = t.test(log1p(.SD$normal.fg[.SD$fused == TRUE]),
                                                      log1p(.SD$normal.fg[.SD$fused == FALSE]))$p.value),
                                         by = .(loose.index)]
            } else {
                ks.cov.dt = data.table(loose.index = numeric(), ks = numeric(), t = numeric())
            }
            
            if ("fused" %in% colnames(mean.cov.dt) & "unfused" %in% colnames(mean.cov.dt)) {
                mean.cov.dt[, norm.change := norm_thresh < abs(fused - unfused)]
            } else {
                mean.cov.dt[, norm.change := NA]
            }

            ## if there is no normal coverage, we cannot confidently tell that there is no change in normal CN
            le.dt$norm.fused = mean.cov.dt$fused[match(le.dt$loose.index, mean.cov.dt$loose.index)]
            le.dt$norm.unfused = mean.cov.dt$unfused[match(le.dt$loose.index, mean.cov.dt$loose.index)]
            le.dt$norm.change = mean.cov.dt$norm.change[match(le.dt$loose.index, mean.cov.dt$loose.index)]

            ## add ks test p values
            le.dt$norm.ks = ks.cov.dt$ks[match(le.dt$loose.index, ks.cov.dt$loose.index)]
            le.dt$norm.t = ks.cov.dt$t[match(le.dt$loose.index, ks.cov.dt$loose.index)]

            ## idk, NA handling?
            le.dt[is.na(norm.fused), norm.change := TRUE]
            le.dt[is.na(norm.unfused), norm.change := TRUE]

            ## track logp value from wilcoxon and t test
            le.dt[, logw := -log10(norm.ks)]
            le.dt[, logt := -log10(norm.t)]
            le.dt[, logfc := log2(norm.fused) - log2(norm.unfused)]
            le.dt[, norm.diff := norm.fused - norm.unfused]
        } else {
            ## if loose ends overlap zero coverage... this is problematic
            le.dt$norm.change = TRUE
            le.dt$norm.fused = NA
            le.dt$norm.unfused = NA
        }
    } else {
        le.dt$norm.change = FALSE
        le.dt$norm.fused = NA
        le.dt$norm.unfused = NA
    }

    ## ## get CN of each side of loose end
    fused.unfused.dt[, fused.cn := gg$nodes$dt$cn[fused.node.id]]
    fused.unfused.dt[, unfused.cn := gg$nodes$dt$cn[unfused.node.id]]
    fused.unfused.dt[, fused.lower := (fused.cn <= unfused.cn)]

    ## merge these with loose ends
    le.dt = merge.data.table(le.dt, fused.unfused.dt, by = "loose.index", all.x = TRUE)

    ## BP mask is still important - the CNP might be absent in this particular matched normal but present
    le.dt[, keep := (!norm.change)]

    ## add stringified loose end
    loose.end.str = gr.string(dt2gr(le.dt[, .(seqnames, start, end, strand)]))
    le.dt[, loose.end := paste0(seqnames, ":", start, strand)]

    return(le.dt)
}

#' @name gg_loose_end
#' @title gg_loose_end
#'
#' @param jabba_rds (character) path to jabba output
#' @param return.type (character) one of GRanges, data.table
#' @param std.chrs (logical) return only loose ends on standard chromosomes
#' 
#' @return loose ends of a JaBbA graph as GRanges
#' if gGraph is not balanced, will produce a warning and return empty GRanges
gg_loose_end = function(jabba_rds = NA_character_,
                        return.type = "GRanges",
                        std.chrs = TRUE) {
    
    if (!check_file(jabba_rds)) {
        stop("file does not exist: ", jabba_rds)
    }
    
    gg = gG(jabba = jabba_rds)

    if (is.null(gg$nodes$dt$loose.cn.left) | is.null(gg$nodes$dt$loose.cn.right)) {
        stop("gg missing field loose.cn.left and loose.cn.right. please check if graph if balanced.")
    }

    ll = gr2dt(gr.start(gg$nodes[!is.na(cn) & loose.cn.left>0]$gr))[, ":="(lcn = loose.cn.left,
                                                                           strand = "+")]
    lr = gr2dt(gr.end(gg$nodes[!is.na(cn) & loose.cn.right>0]$gr))[, ":="(lcn = loose.cn.right,
                                                                           strand = "-")]
    if (nrow(ll)) {
        ll[, ":="(lcn = loose.cn.left, strand = "+")]
    }

    if (nrow(lr)) {
        lr[, ":="(lcn = loose.cn.right, strand = "-")]
    }

    all.loose.dt = rbind(ll, lr, fill = TRUE)

    if (all.loose.dt[, .N]) {
        stdchrs = c(as.character(1:22), "X", "Y")
        all.loose.dt = all.loose.dt[(as.character(seqnames) %in% stdchrs) |
                                    (as.character(seqnames) %in% paste0("chr", stdchrs)),]
    }


    if (return.type == "GRanges") {
        all.loose.gr = dt2gr(all.loose.dt, seqlengths = seqlengths(gg))
        return(all.loose.gr)
    }
    return(all.loose.dt)
}
