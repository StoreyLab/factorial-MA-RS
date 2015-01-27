#' create R datasets that end up in the data/ subdirectory,
#' such as the RNA-Seq count matrix and design.
#'
#' @details Note that microarrayMatrix.rda is not created by this
#' function (has to be pre-included)
#'
#' @param outfolder Folder to write .rda files to
#'
#' @importFrom tidyr separate
#'
#' @examples
#'
#' \dontrun{
#' # run from inside the factorial folder to regenerate R data
#' preprocess("data")
#' }
#'
#' @export
preprocess <- function(outfolder = "data") {
    dir.create(outfolder, showWarnings = FALSE)

    # load RNA-Seq matrix
    data_file <- function(f) system.file("extdata", f, package = "factorial")
    out_file <- function(n) paste(outfolder, paste0(n, ".rda"), sep="/")

    raw.counts <- read.table(data_file("count_matrix.txt"), sep = "\t",
                                 header = TRUE, row.names = 1)

    experimentDesign <- do.call(rbind, strsplit(colnames(raw.counts), ""))
    experimentDesign <- as.data.frame(experimentDesign)
    colnames(experimentDesign) <- c("condition", "biological",
                                    "preparation", "lane")

    # reorder design and raw counts
    ord <- do.call(order, experimentDesign)
    experimentDesign <- experimentDesign[ord, ]
    raw.counts <- raw.counts[ord]

    rnaseqMatrix <- as.matrix(raw.counts[grep("^Y", rownames(raw.counts)), ])
    statistic.counts <- raw.counts[grep("^[A-Z]", rownames(raw.counts), invert = TRUE), ]

    # save to outfolder
    save(rnaseqMatrix, file = out_file("rnaseqMatrix"))
    save(experimentDesign, file = out_file("experimentDesign"))

    # turn microarray.Rdata (created somewhere else) into appropriate objects
    load(data_file("microarray.Rdata"))

    microarrayRedMatrix <- as.matrix(red.matrix)
    microarrayGreenMatrix <- as.matrix(green.matrix)
    microarrayMatrix <- log2(microarrayRedMatrix / microarrayGreenMatrix)
    save(microarrayRedMatrix, file = out_file("microarrayRedMatrix"))
    save(microarrayGreenMatrix, file = out_file("microarrayGreenMatrix"))
    save(microarrayMatrix, file = out_file("microarrayMatrix"))

    # preprocess qPCR data
    qPCR.raw <- read.csv(data_file("qPCR_data.txt"), sep = "\t",
                         na = c("Undetermined", "neg control"))

    qPCR.raw <- qPCR.raw %>% mutate(technical=factor(technical))

    # one run of rRNA gave problematic results. We re-ran it and it worked: the
    # problematic results are contained as rRNA18S-X
    rRNA.rejected <- qPCR.raw %>% filter(gene == "rRNA18S-X")

    # not used, so not saved
    rRNA <- qPCR.raw %>% filter(gene == "rRNA18S") %>% filter(!is.na(ct))
    qPCR <- qPCR.raw %>% filter(!grepl("rRNA18S", gene))

    qPCR <- qPCR %>% tidyr::separate(sample, c("condition", "biological"), 1) %>%
        filter(!is.na(ct)) %>% dplyr::select(-dilution) %>%
        rename(symbol = gene) %>%
        mutate(gene = as.character(revmap(org.Sc.sgdGENENAME))[as.character(symbol)])
    save(qPCR, file = out_file("qPCR"))

    # save platereader concentrations
    plateReader <- read.delim(data_file("plate_reader_09092014.txt"))
    save(plateReader, file = out_file("plateReader"))

    # save dubious genes
    description = as.character(org.Sc.sgd.db::org.Sc.sgdDESCRIPTION)
    dubious = names(description)[grepl("Dubious", description)]
    save(dubious, file = out_file("dubious"))

    # save pathway genes
    pathway <- read.delim(data_file("special_genes.txt"))
    save(pathway, file = out_file("pathway"))

    # some extra stuff:
    # SRA descriptions:
    # cat(with(experimentDesign, paste0(colnames(rnaseqMatrix), " - ", "Condition: ", c("ethanol", "glucose")[condition], ", Biological Replicate: ", biological, ", Preparation: ", preparation, ", ", "Lane: ", lane)), sep = "\n")
}
