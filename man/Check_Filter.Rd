\name{Check_Filter}
\alias{Check_Filter}
\title{Check different filtering cut-offs for RNAseq data.}
\usage{
Check_Filter(D, thresh=0)
}
\description{
Plots a density plot of summed raw counts per gene at various filtering thresholds.}
\arguments{
\item{D}{A \code{DESeqDataSet} object}
\item{thresh}{numeric, filtering threshold to test}
}
\value{
Returns a plot.
}
\examples{
# Check_Filter(D, thresh=0)
}
