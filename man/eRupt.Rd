\name{eRupt}
\alias{eRupt}
\title{Makes a volcano plot.}
\usage{
eRupt(results, alpha = 0.05, FCsig = 1)
}
\description{
Make a volcano plot from a \code{makeComp()} output file.}
\arguments{
\item{results}{\code{data.frame}, a \code{makeComp()} output file}
\item{alpha}{\code{numeric}, pAdj threshold for significance}
\item{FCsig}{\code{numeric}, absolute log2FC threshold for significance}
}
\value{
Returns a volcano plot.
}
\examples{
# eRupt(results)
}
