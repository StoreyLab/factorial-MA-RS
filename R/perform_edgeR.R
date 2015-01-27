#' Perform edgeR exact test on a matrix, given the groups
#'
#' @import edgeR
#' @import biobroom
#'
#' @export
perform_edgeR <- function(matrix, group=c("E", "E", "G", "G")) {
    dge.list <- DGEList(counts=matrix, group=group)
    dge.list <- calcNormFactors(dge.list)
    dge.list <- estimateCommonDisp(dge.list)
    dge.list <- estimateTagwiseDisp(dge.list)
    ret <- exactTest(dge.list, pair=c("E", "G"))
    ret <- biobroom::tidy(ret)

    # add additional interesting columns
    ret$intensity <- rowSums(matrix)
    ret$intensityE <- rowSums(matrix[, group == "E"])
    ret$intensityG <- rowSums(matrix[, group == "G"])
    ret
}
