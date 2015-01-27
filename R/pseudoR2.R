
#' @export
design_matrices <- function() {
    mm = model.matrix(~ condition + condition:biological + condition:biological:preparation, experimentDesign)

    dfs = lapply(c(1, 2, 4, 8), function(k) 1:k)
    design.matrices = c(lapply(dfs, function(k) mm[, k]))
    design.matrices[[1]] = matrix(design.matrices[[1]])
    design.matrices
}

# n = NCOL(gene.counts.filtered)
# p = c(1, 2, 4, 8)

#' calculate pseudo-R^2 for every row based on a dispersion model
#'
#' @param mat count matrix
#' @param disp Calculated per-gene dispersion, or 0 for a Poisson model
#'
#' @import data.table
#' @import edgeR
#' @importFrom reshape2 melt
#'
#' @details This function is largely from the pre-refactoring version and could
#' use some further reorganization.
#'
#' @export
constructPseudoR2 = function(mat, disp) {
    #fit = glmFit(d, mm, dispersion=disp)
    mm = model.matrix(~ condition + condition:biological + condition:biological:preparation, experimentDesign)
    d <- DGEList(counts=mat, group=experimentDesign$condition)
    d <- calcNormFactors(d)
    d <- estimateGLMTrendedDisp(d, mm)
    d <- estimateGLMTagwiseDisp(d, mm)

    fits = sapply(design_matrices(), function(m) edgeR::glmFit(d, m, dispersion=disp))
    fitted.vals.each = lapply(fits, fitted.values)

    y = mat
    ybar = fitted.vals.each[[1]]

    p = c(NA, 1, 3, 7)
    n <- ncol(mat)

    disp.inv = 1 / disp
    R2s.full = sapply(2:length(fitted.vals.each), function(i) {
        mu = fitted.vals.each[[i]]
        if (length(disp) == 1 & disp[1] == 0) {
            num = rowSums(y * log(y / mu) - (y - mu), na.rm=TRUE)
            denom = rowSums(y * log(y / ybar) - (y - ybar), na.rm=TRUE)
        } else {
            num = rowSums(y * log(y / mu) - (y + disp.inv) * log((y + disp.inv) / (mu + disp.inv)), na.rm=TRUE)
            denom = rowSums(y * log(y / ybar) - (y + disp.inv) * log((y + disp.inv) / (ybar + disp.inv)), na.rm=TRUE)
        }
        (1 - (num * (n - 1)) / (denom * (n - p[i] - 1)))
    })
    # adjust
    #p = c(1, 3, 7)
    #R2s.adj = t(1 - (1 - t(R2s.full) * (n - 1) / (n - p - 1)))
    adj.pseudo.R2s = t(apply(R2s.full, 1, function(r) diff(c(0, r, 1))))
    #logliks.each = sapply(fitted.vals.each, function(m) logliks(counts.m, m, disp))
    #diffs = t(apply(logliks.each, 1, diff))
    #pseudo.R2s = (diffs / rowSums(diffs))
    #adj.pseudo.R2s = t(1 - (1 - t(pseudo.R2s) * (n - 1) / (n - p - 1)))

    #adj.pseudo.R2s = t(1 - (1 - t(pseudo.R2s) * (n - 1) / (n - p - 1)))
    adj.pseudo.R2s = as.data.frame(adj.pseudo.R2s)
    colnames(adj.pseudo.R2s) = colnames(experimentDesign)
    adj.pseudo.R2s$Depth = rowSums(mat)
    adj.pseudo.R2s$Gene = rownames(mat)

    # remove the NA cases
    adj.pseudo.R2s = adj.pseudo.R2s[complete.cases(adj.pseudo.R2s), ]

    pseudo.R2s.m = as.data.table(melt(adj.pseudo.R2s, id=c("Gene", "Depth")))
    setnames(pseudo.R2s.m, "variable", "Level")
    pseudo.R2s.m[, Quantile:=rank(Depth) / length(Depth), by="Level"]

    #list(fitted.vals.each=fitted.vals.each, logliks.each=logliks.each,
    #     pseudo.R2s=pseudo.R2s, adj.pseudo.R2s=adj.pseudo.R2s,
    #     pseudo.R2s.m=pseudo.R2s.m)
    pseudo.R2s.m
}
