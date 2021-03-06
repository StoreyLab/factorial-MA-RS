% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/preprocess.R
\name{preprocess}
\alias{preprocess}
\title{create R datasets that end up in the data/ subdirectory,
such as the RNA-Seq count matrix and design.}
\usage{
preprocess(outfolder = "data")
}
\arguments{
\item{outfolder}{Folder to write .rda files to}
}
\description{
create R datasets that end up in the data/ subdirectory,
such as the RNA-Seq count matrix and design.
}
\details{
Note that microarrayMatrix.rda is not created by this
function (has to be pre-included)
}
\examples{
\dontrun{
# run from inside the factorial folder to regenerate R data
preprocess("data")
}
}

