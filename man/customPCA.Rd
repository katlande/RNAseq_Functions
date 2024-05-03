\name{customPCA}
\alias{customPCA}
\title{Plots a PCA of samples.}
\usage{
customPCA(D, meta, ptcol, filter_threshold=5,  ptsize=2,
         labsize=3, ptshape=NA, main="PCA", label=F, facet=NA,
         frows=NA, fcols=NA, fscale="fixed", custom_cols=NULL)
}
\description{
Checks if meta data and count data match; reorders meta data when necessary, and warns users if files do not match.}
\arguments{
\item{D}{A \code{DESeqDataSet} object}
\item{meta}{A \code{data.frame} of meta data with sample names as \code{row.names}}
\item{ptcol}{\code{character}, name of variable in \code{meta} to define point colour}
\item{filter_threshold}{\code{numeric}, minimum counts/gene for gene inclusion}
\item{ptsize}{\code{numeric}, size of plotted points}
\item{labsize}{\code{numeric}, size of label text}
\item{ptshape}{\code{character}, name of variable in \code{meta} to define point shape, if not \code{NA}}
\item{main}{plot title}
\item{label}{\code{boolean}}
\item{facet}{\code{character}, name of variable in \code{meta} to define plot facets, if not \code{NA}}
\item{frows}{if \code{facet} is not \code{NA}, number of facet rows}
\item{fcols}{if \code{facet} is not \code{NA}, number of facet columns}
\item{fscale}{if \code{facet} is not \code{NA}, facet scales: "fixed", "free", "free_x", "free_y"}
\item{custom_cols}{\code{vector} of colours for point colours}
}
\value{
Returns PCA plot of data.
}
\examples{
# customPCA(D, meta, ptcol="Age", ptshape="Sex")
}
