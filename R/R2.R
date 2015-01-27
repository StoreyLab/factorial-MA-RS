#' compute R-squared values for many genes
#'
#' @param m a normalized matrix from normalizeBetweenArrays or voom
#' @param mm.full the model matrix whose R^2 is being examined
#' @param adjusted whether the R^2 values should be adjusted for degrees of freedom
#'
#' @export
R2s = function(m, mm.full, adjusted=TRUE) {
    # given a matrix and model matrices for full and null
    lm.full = limma::lmFit(m, mm.full)
    RSS.full = rowSums(residuals(lm.full, m)^2)
    SS.total = rowSums((m - rowMeans(m))^2)
    n = NCOL(m)
    if (adjusted) {
        (1 - ((RSS.full / lm.full$df.residual) / (SS.total / (n - 1))))
    }
    else {
        (1 - (RSS.full / SS.total))
    }
}

#' construct R2 across all conditions
#'
#' @param m normalized matrix
#' @param intensity vector of intensities
#'
#' @importFrom plyr revalue
#' @import dplyr
#'
#' @export
construct_R2s = function(m, intensity) {
    data(experimentDesign)
    colnames(experimentDesign)[4] = "chip.lane"

    mms <- with(experimentDesign, list(model.matrix(~ condition),#,
                                      model.matrix(~ 0 + condition:biological),
                                      model.matrix(~ 0 + condition:biological:preparation)))
    ret <- do.call(cbind, lapply(mms, function(mm) R2s(m, mm)))
    # take cumulative differences, add residuals
    ret <- t(apply(ret, 1, function(r) diff(c(0, r, 1))))

    colnames(ret) <- c("Condition", "Biological", "Preparation", "Chip.Lane")
    ret <- data.frame(cbind(ret, Intensity=intensity))
    ret$Gene <- rownames(m)

    ret <- ret %>% tidyr::gather(Level, value, -Intensity, -Gene)
    ret$Level <- plyr::revalue(ret$Level, c(Chip.Lane="Chip/Lane"))
    ret
}
