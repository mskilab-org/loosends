dflm = function (x, last = FALSE, nm = "") 
{
    if (is.null(x)) 
        out = data.frame(name = nm, method = as.character(NA), 
            p = as.numeric(NA), estimate = as.numeric(NA), ci.lower = as.numeric(NA), 
            ci.upper = as.numeric(NA), effect = as.character(NA))
    else if (any(c("lm", "betareg") %in% class(x))) {
        coef = as.data.frame(summary(x)$coefficients)
        colnames(coef) = c("estimate", "se", "stat", "p")
        if (last) 
            coef = coef[nrow(coef), ]
        coef$ci.lower = coef$estimate - 1.96 * coef$se
        coef$ci.upper = coef$estimate + 1.96 * coef$se
        if (!is.null(summary(x)$family)) {
            fam = summary(x)$family$family
            if (summary(x)$family$link %in% c("log", "logit")) {
                coef$estimate = exp(coef$estimate)
                coef$ci.upper = exp(coef$ci.upper)
                coef$ci.lower = exp(coef$ci.lower)
            }
        }
        else fam = "Unknown"
        if (!last) 
            nm = paste(nm, rownames(coef))
        out = data.frame(name = nm, method = fam, p = signif(coef$p, 
            3), estimate = coef$estimate, ci.lower = coef$ci.lower, 
            ci.upper = coef$ci.upper, effect = paste(signif(coef$estimate, 
                3), " [", signif(coef$ci.lower, 3), "-", signif(coef$ci.upper, 
                3), "]", sep = ""))
    }
    else if (class(x) == "htest") {
        if (is.null(x$estimate)) 
            x$estimate = x$statistic
        if (is.null(x$conf.int)) 
            x$conf.int = c(NA, NA)
        out = data.table(name = nm, method = x$method, estimate = x$estimate, 
            ci.lower = x$conf.int[1], ci.upper = x$conf.int[2], 
            effect = paste(signif(x$estimate, 3), " [", signif(x$conf.int[1], 
                3), "-", signif(x$conf.int[2], 3), "]", sep = ""), 
            p = x$p.value)
    }
    else if (class(x) == "polr") {
        coef = coef(summary(x)) %>% as.data.frame
        nm = paste(nm, rownames(coef))
        coef = as.data.table(coef)
        setnames(coef, c("estimate", "se", "t"))
        out = data.table(name = nm) %>% cbind(coef)
        out$p = pnorm(abs(out$t), lower.tail = FALSE) * 2
        out[, `:=`(ci.lower, estimate - 1.96 * se)]
        out[, `:=`(ci.upper, estimate + 1.96 * se)]
        out[, `:=`(effect, paste(signif(estimate, 3), " [", signif(ci.lower, 
            3), "-", signif(ci.upper, 3), "]", sep = ""))]
    }
    else {
        out = data.frame(name = nm, method = x$method, p = signif(x$p.value, 
            3), estimate = x$estimate, ci.lower = x$conf.int[1], 
            ci.upper = x$conf.int[2], effect = paste(signif(x$estimate, 
                3), " [", signif(x$conf.int[1], 3), "-", signif(x$conf.int[2], 
                3), "]", sep = ""))
    }
    out$effect = as.character(out$effect)
    out$name = as.character(out$name)
    out$method = as.character(out$method)
    rownames(out) = NULL
    return(as.data.table(out))
}

is.dup = function (x) 
{
    if (is.matrix(x)) 
        x = as.data.frame(x)
    if (is.data.frame(x) | data.table::is.data.table(x)) {
        tmp = x[[1]]
        if (ncol(x) > 1) 
            for (i in 2:ncol(x)) tmp = paste(tmp, x[[i]], sep = "@!$!$!@")
        x = tmp
    }
    d = duplicated(x)
    return(x %in% x[d])
}
