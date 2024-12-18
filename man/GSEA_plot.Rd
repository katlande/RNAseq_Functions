\name{GSEA_plot}
\alias{GSEA_plot}
\title{Plot GSEA-mode WebGestaltPipe output files.}
\usage{
GSEA_plot(GSEA, balance=T, n=10, ontology_name="KEGG", logSize=F)
}
\description{
Plots a dotplot of top ORA hits.}
\arguments{
\item{GSEA}{\code{data.frame}, WebGestaltPipe() output file.}
\item{balance}{boolean, whether or not to include an equal number of up and down regulated terms.}
\item{n}{numeric, total number of terms to use.}
\item{ontology_name}{character, ontology name for plot output.}
\item{logSize}{boolean, whether to log10-scale the point size parameter.}
}
\value{
Returns a plot.
}
\examples{
# GSEA_plot(GSEA)
}
