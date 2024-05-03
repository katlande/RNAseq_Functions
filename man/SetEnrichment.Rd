\name{SetEnrichment}
\alias{SetEnrichment}
\title{Check if a gene set is enriched in another gene set.}
\usage{
SetEnrichment(geneList, bkgd_subset, bkgd)
}
\description{
Pipe to a fisher's exact test that checks if genes in geneList appear in bkgd_subset more often than expected by random chance.}
\arguments{
\item{geneList}{\code{character}, gene set of interest}
\item{bkgd_subset}{\code{character}, gene set to query against (e.g., differentially expressed genes)}
\item{bkgd}{\code{character}, background (e.g., all filtered genes))}
}
\value{
Returns a \code{vector} containing the p-value of enrichment at index 1, and the enrichment ratio at index 2.
}
\examples{
# SetEnrichment(GOIs, DEGs, allGenes)
}
