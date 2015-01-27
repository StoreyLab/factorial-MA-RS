#' Given a microarray/RNASeq matrix and a subset of columns, pool within columns
#'
#' @param m a matrix of RNA-Seq counts or microarray logFCs
#' @param des Experiment design (such as built-in experimentDesign)
#' @param cond vector used to subset expression matrix and design
#' @param average Whether values should be averaged: otherwise they are summed
#'
#' @return The matrix, using the given subset, pooled within condition:bio
#'
#' @export
pool_matrix.matrix <- function(m, des, cond, average = FALSE) {
    d <- des[cond, ]
    mm <- model.matrix(~ biological:condition, d)[, -1]
    sm <- as.matrix(m[, cond])

    pooled <- sm %*% mm
    if (average) {
        pooled <- t(t(pooled) / colSums(mm))
    }

    n <- colnames(mm)
    colnames(pooled) <- paste0(substr(n, nchar(n), nchar(n)), 1:2)
    pooled
}

#' method for pooling an RGList object
#'
#' @export
pool_matrix.RGList <- function(m, des, cond, average = TRUE) {
    # pool R and G separately, each as geometric means
    poolR <- exp(pool_matrix(log(m$R), des, cond, average))
    poolG <- exp(pool_matrix(log(m$G), des, cond, average))

    new("RGList", list(R=poolR, G=poolG))
}

#' @export
pool_matrix <- function(m, des, cond, average = FALSE) UseMethod("pool_matrix")
