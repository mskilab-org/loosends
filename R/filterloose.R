#' @name filter_loose
#' @title filter_loose
#'
#' @param jabba_rds (character) path to JaBbA
#' @param coverage (character) path to coverage file with gc-adjusted counts
#' @param background (character) path to dryclean with background
#' @param bg.field (character) background coverage field (background)
#' @param norm.field (character) normal count field (reads.corrected)
#' @param mask (character) path to coverage mask
#' @param bp.mask (character) path to coverage mask
#' @param mask_pad (numeric) pad on distance to mask (in practice just for annotation purposes)
#' @param norm_thresh (numeric) 0.6 (from filter.loose)
#' @param bg_thresh (numeric) 0.9 threshold for germline CN change
#' @param min_size (numeric) minimum size in BP
#' @param subset (logical) subset loose ends to include only valid ones before returning
#' @param verbose (logical)
grab_loose = function(jabba_rds,
                      coverage = "/dev/null",
                      background = "/dev/null",
                      bg.field = "background",
                      norm.field = 'reads.corrected',
                      ## cov.pad = 5e4,
                      mask = "~/projects/gGnome/files/zc_stash/tile.mask.rds",
                      bp.mask = "/dev/null", ##"~/projects/gGnome/files/zc_stash/common.bps.rds",
                      mask_pad = 0,
                      norm_thresh = 0.6,
                      bg_thresh = 0.9,
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

    if (!check_file(background)) {
        if (verbose) {
            message("Background file not supplied!")
        }
        use.bg = FALSE
    } else {
        if (verbose) {
            message("Reading background coverage")
        }
        bg.gr = readRDS(background)
        bg.gr = gr.nochr(bg.gr)
        use.bg = TRUE
    }

    if (use.coverage) {
        if (!(norm.field %in% names(values(cov.gr)))) {
            stop("supplied field not in coverage metadata: ", norm.field)
        }
        cov.gr = cov.gr %Q% (!is.infinite(values(cov.gr)[[norm.field]]))
    }

    if (use.bg) {
        if (!(bg.field %in% names(values(bg.gr)))) {
            stop("supplied field not in coverage metadata: ", bg.field)
        }
        bg.gr = bg.gr %Q% (!is.infinite(values(bg.gr)[[bg.field]]))
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

    if (check_file(bp.mask)) {

        if (verbose) {
            message("Checking overlaps with common breakpoints!")
        }
        
        bp.mask.gr = readRDS(bp.mask)
        le.dt[, bp.mask := le.gr %^% gr.stripstrand(bp.mask.gr)]
    } else {
        le.dt[, bp.mask := FALSE]
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
        cov.gr$abscn = skitools::rel2abs(cov.gr, field = norm.field, purity = 1, ploidy = 2)
    }

    if (use.bg) {
        if (verbose) { message("rel2abs transforming background coverage") }
        ##mean.ploidy = (gg$meta$purity * gg$meta$ploidy) + (1 - gg$meta$purity) * 2
        bg.gr$abscn = skitools::rel2abs(bg.gr, field = bg.field, purity = 1, ploidy = 2)
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

    if (use.bg) {

        ## add background segmentation here, which is hopefully not too noisy
        bg.cov.dt = gr.findoverlaps(bg.gr, segs.gr, return.type = 'data.table')
        if (nrow(sides.cov.dt)) {
            bg.cov.dt[, bg.cov := values(bg.gr)[["abscn"]][query.id]]
            bg.cov.dt[, fused := values(segs.gr)$fused[subject.id]]
            bg.cov.dt[, loose.index := values(segs.gr)$loose.index[subject.id]]
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
                                         .(## fused.sw = shapiro.test(log1p(.SD$normal.fg[.SD$fused == TRUE]))$p.value,
                                           ## unfused.sw = shapiro.test(log1p(.SD$normal.fg[.SD$fused == FALSE]))$p.value,
                                           ks = wilcox.test(.SD$normal.fg[.SD$fused == TRUE],
                                                            .SD$normal.fg[.SD$fused == FALSE])$p.value,
                                           t = t.test(log1p(.SD$normal.fg[.SD$fused == TRUE]),
                                                      log1p(.SD$normal.fg[.SD$fused == FALSE]))$p.value),
                                         by = .(loose.index)]
            } else {
                ks.cov.dt = data.table(loose.index = numeric(), ks = numeric(), t = numeric())
            }
            
            if ("fused" %in% colnames(mean.cov.dt) & "unfused" %in% colnames(mean.cov.dt)) {
                ## should be rel2abs transformed
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
            ## le.dt$fused.sw = ks.cov.dt$fused.sw[match(le.dt$loose.index, ks.cov.dt$loose.index)]
            ## le.dt$unfused.sw = ks.cov.dt$unfused.sw[match(le.dt$loose.index, ks.cov.dt$loose.index)]

            ## idk, NA handling?
            le.dt[is.na(norm.fused), norm.change := TRUE]
            le.dt[is.na(norm.unfused), norm.change := TRUE]

            ## track logp value from wilcoxon and t test
            ## track log fold change and difference between fused and unfused
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

    if (use.bg) {

        ## compute mean coverage for fused and unfused sides
        if (sides.cov.dt[, .N]) {
            mean.cov.dt = bg.cov.dt[, .(mean.normal.cov = mean(bg.cov, na.rm = T)),
                                       by = .(loose.index, fused)]
            mean.cov.dt[, fused := ifelse(fused, 'fused', 'unfused')]
            mean.cov.dt[, fused := factor(fused, levels = c("fused", "unfused"))]
            mean.cov.dt = dcast.data.table(mean.cov.dt, loose.index ~ fused, value.var = "mean.normal.cov")
            
            ## determine whether it's higher than normal beta
            if ("fused" %in% colnames(mean.cov.dt) & "unfused" %in% colnames(mean.cov.dt)) {
                mean.cov.dt[, bg.change := bg_thresh < abs(fused - unfused)]
            } else {
                mean.cov.dt[, bg.change := NA]
            }

            ## mark normal change and copy normal fused and unfused coverage
            le.dt$bg.change = le.dt$loose.index %in% mean.cov.dt[(bg.change), loose.index]

            ## if there is no normal coverage, we cannot confidently tell that there is no change in normal CN
            le.dt$bg.fused = mean.cov.dt$fused[match(le.dt$loose.index, mean.cov.dt$loose.index)]
            le.dt$bg.unfused = mean.cov.dt$unfused[match(le.dt$loose.index, mean.cov.dt$loose.index)]
            le.dt[is.na(bg.fused), bg.change := TRUE]
            le.dt[is.na(bg.unfused), bg.change := TRUE]
        } else {
            le.dt$bg.change = TRUE
            le.dt$bg.fused = NA
            le.dt$bg.unfused = NA
        }
    } else {
        le.dt$bg.change = FALSE
        le.dt$bg.fused = NA
        le.dt$bg.unfused = NA
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
