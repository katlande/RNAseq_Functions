\name{Reliable_heatmap}
\alias{Reliable_heatmap}
\title{A pipe to a reproducible aesthetic for pheatmap.}
\usage{
Reliable_heatmap(norm, genes=NULL, sample_meta=NULL, annotation=F, annots=NULL, title="heatmap", 
nrow=1, ncol=1, log=F, angle="45", frow=8, fcol=8, annots_as_factors=T, rowlabs=T, collabs=T,
cols=c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")
}
\description{Makes a pheatmap with a reproducible aesthetic.}
\arguments{
\item{norm}{\code{matrix}, numeric matrix of input data}
\item{genes}{\code{character}, optional subsetting vector for \code{norm}; will only plot rows with \code{row.names} found in \code{genes}}
\item{sample_meta}{\code{data.frame}, optional meta data for norm, where \code{row.names} of sample_meta match \code{colnames} of norm. Function will only plot samples found in \code{sample_meta}, thus this option can be used as a way to subset which columns of norm are plotted. \code{sample_meta} must be supplied if annotation=T}
\item{annotation}{\code{boolean}, whether to annotate heatmap columns by \code{sample_meta} columns}
\item{annots}{\code{character}, optional vector of \code{sample_meta} colnames to use for annotation}
\item{title}{\code{character}, heatmap title}
\item{nrow}{\code{numeric}, pipe to \code{pheatmap} cuttree_rows}
\item{ncol}{\code{numeric}, pipe to \code{pheatmap} cuttree_cols}
\item{log}{\code{boolean}, whether to take the natural log of the input data}
\item{angle}{\code{character}, pipe to \code{pheatmap} angle_col: “0”, “45”, “90”, “270”, “315” }
\item{frow}{\code{numeric}, fontsize of row labels}
\item{fcol}{\code{numeric}, fontsize of column labels}
\item{annots_as_factors}{\code{boolean}, if true, numeric annotation columns will be labeled as factors}
\item{rowlabs}{\code{boolean}, whether to plot row labels}
\item{collabs}{\code{boolean}, whether to plot column labels}
\item{cols}{\code{character}, custom vector of colours to use for cell colours}
\value{
Returns a pheatmap object.
}
\examples{
# Reliable_heatmap(norm)
# Reliable_heatmap(norm, genes=c("gene1", "gene2", ... "geneN"), sample_meta=meta_data, annotation=T)
}

