\name{CheckDepth}
\alias{CheckDepth}
\title{Checks the read depth of all samples.}
\usage{
CheckDepth(raw, minDepth=1e07)
}
\description{
Plots boxplot of raw read depth from all samples.}
\arguments{
\item{raw}{A gene x sample \code{data.frame} of raw reads}
\item{minDepth}{numeric, minimum depth cut-off (1e07 is recommended transcriptomes similar in size to mouse, human).}
}
\value{
Returns a plot.
}
\examples{
# CheckDepth(raw)
}
