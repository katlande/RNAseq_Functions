\name{Enrich2Table}
\alias{Enrich2Table}
\title{Convert top WebGestalt_Pipe results to a table grob.}
\usage{
Enrich2Table(df, mode="ORA", n=10)
}
\description{
Plots a table with enrichment result data.}
\arguments{
\item{df}{A \code{data.frame} output by WebGestalt_Pipe()}
\item{mode}{"ORA" or "GSEA", depending on which mode was used in WebGestalt_Pipe(), Auto="ORA"}
\item{n}{\code{integer}, number of rows to print in output table, Auto=10}
}
\value{
Returns a plot.
}
\examples{
# Check_Filter(D, thresh=0)
}
