\name{ORA_plot}
\alias{ORA_plot}
\title{Plot ORA-mode WebGestaltPipe output files.}
\usage{
ORA_plot(ORA, balance=T, n=10, ontology_name="Biological Process", logSize=F)
}
\description{
Plots a dotplot of top ORA hits.}
\arguments{
\item{ORA}{\code{data.frame}, WebGestaltPipe() output file.}
\item{balance}{boolean, whether or not to include an equal number of up and down regulated terms.}
\item{n}{numeric, total number of terms to use.}
\item{ontology_name}{character, ontology name for pot output.}
\item{logSize}{boolean, whether to log10-scale the point size parameter.}
}
\value{
Returns a plot.
}
\examples{
# ORA_plot(ORA)
}
