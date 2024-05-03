DESeq_and_Filt <- function(D, filter_thresh, minSamples=3){
  keep <- rowSums(counts(D)) >= filter_thresh & rowSums(counts(D) > 0) > minSamples # filter
  D <- D[keep,]
  D <- DESeq2::DESeq(D) # find DEGs
  return(D)
}
