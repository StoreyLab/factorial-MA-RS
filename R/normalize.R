#' Normalize an RS/MA matrix according to the appropriate normalization method
#'
#' @export
normalizeMatrix <- function(mat, normalize.method="scale") UseMethod("normalizeMatrix")

#' @import limma
#'
#' @export
normalizeMatrix.RGList <- function(mat, normalize.method="scale") {
    if (normalize.method %in% c("TMM", "RLE")) {
        return(NULL)
    }

    limma::normalizeBetweenArrays(mat, method=normalize.method)
}


#' Normalize an RNA-Seq matrix
#'
#' @export
normalizeMatrix.matrix <- function(mat, normalize.method="scale") {
    if (normalize.method %in% c("TMM", "RLE")) {
        TMM.factors = colSums(mat) * edgeR::calcNormFactors(mat, method=normalize.method)
        mat = log(t(t(mat + 1) / TMM.factors))
    } else {
        mat <- limma::voom(mat, normalize.method=normalize.method)
    }
    mat
}
