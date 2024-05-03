\name{SetScores}
\alias{SetScores}
\title{Get the median Z-score of normalized counts for a set of genes in all samples}
\usage{
SetScores(norm, goi, meta)
}
\description{Scales the normalized counts matrix across rows and takes the median value of all genes in a specified set of genes for each sample.}
\arguments{
\item{norm}{\code{data.frame}, gene x sample matrix of normalized counts}
\item{goi}{\code{character} or \code{list}. If querying one gene set, provide it as a \code{character} vector. If querying multiple gene sets, provide them as a named \code{list} of \code{character} vectors.}
\item{meta}{\code{data.frame} of meta data for all samples in \code{norm}}
}
\value{
Returns a \code{data.frame} containing scores for each sample and input gene set.
}
\examples{
# SetScores(norm, c("Stat1", "Myc", "Ccnd1"), meta)
}
