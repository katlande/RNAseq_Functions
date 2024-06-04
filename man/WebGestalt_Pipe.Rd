\name{WebGestalt_Pipe}
\alias{WebGestalt_Pipe}
\title{Run enrichment analysis through WebGestalt.}
\usage{
WebGestalt_Pipe(mode="ORA", results, alpha=0.05, FC=0.5, DB="geneontology_Biological_Process_noRedundant", species="mmusculus", projectName="scratch")
}
\description{
A pipe that runs a simple ORA or GSEA on a makeComp() output file.}
\arguments{
\item{mode}{"ORA" or "GSEA". Note: GSEA mode runs GSEA on the log2FC values.}
\item{results}{A \code{makeComp()} output \code{data.frame}. Any data frame with Gene name, pAdj, and log2FC can also work. Gene names must be in column 7, and log2foldchange must be in column 2. Columns must be named: "Gene", "padj", and "log2FoldChange".}
\item{alpha}{pAdj cut-off for DEGs (ORA only)}
\item{FC}{absolute log2FC cut-off for DEGs (ORA only)}
\item{DB}{WebGestalt enrichDatabase to query against, see \code{WebGestaltR::listGeneSet()} for a full list.}
\item{species}{WebGestalt organism, see \code{WebGestaltR::listOrganism()} for a full list.}
\item{projectName}{Name of raw data folder. If NULL, no folder is made.}
}
\value{
Returns a \code{data.frame} containing enrichment results.
}
\examples{
# ORA <- WebGestalt_Pipe("ORA", results)
}
