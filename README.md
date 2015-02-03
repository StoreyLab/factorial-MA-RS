A nested parallel experiment demonstrates differences in intensity-dependence between RNA-seq and microarrays
====================

Draft of the manuscript [A nested parallel experiment demonstrates differences in intensity-dependence between RNA-seq and microarrays](http://biorxiv.org/content/early/2014/12/30/013342), developed as a package (see [here](http://rmflight.github.io/posts/2014/07/vignetteAnalysis.html)). Note:

* Draft source is in vignettes, compiled draft is in `inst/doc`
    * Makes use of functions in `R/`
* Raw data is in `inst/extdata`, relevant R data files are in `data`
    * The `preprocess` function, when run within the main directory, turns the data in `inst/extdata` into the `.rda` files in `data`. The vignette (manuscript) makes use only of the `.rda` files built into the package
    * Some datasets have documentation, others don't yet but will

Installation
--------------

You'll first need to install the following packages:

    source("http://bioconductor.org/biocLite.R")
    biocLite(c("limma", "edgeR", "DESeq2", "org.Sc.sgd.db", "GO.db"))
    
    install.packages("devtools")
    devtools::install_github(c("dgrtwo/broom", "dgrtwo/biobroom", "dgrtwo/GSEAMA"))

Then open `factorial.Rproj` and under the Build tab, click "Build & Reload." It should install the remaining packages. You can then go to `vignettes/factorial.Rnw` and click Compile PDF.

To compile the supplemental material, go to `vignettes/supplemental.Rnw` and click Compile PDF again.

Deposited Data
--------------

The microarray data used in this experiment is available from the Gene Expression Omnibus (accession [GSE65365](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65365)), and the RNA-Seq reads are available from the NCBI Short Read Archive (BioProject accession [PRJNA271248](http://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA271248)). 
