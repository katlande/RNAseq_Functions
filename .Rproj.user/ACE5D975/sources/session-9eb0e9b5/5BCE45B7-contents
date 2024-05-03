filt.norm <- function(D, thresh=5, return.tab=T){
  keep <- rowSums(counts(D)) >= thresh # filter raw counts
  D <- D[keep,]
  D <- DESeq2::estimateSizeFactors(D) # normalize data

  if(return.tab==T){
    normalized_counts <- as.data.frame(counts(D, normalized=TRUE))
    return(normalized_counts) # return data as a data frame
  }else{
    return(D)  # return data as a  DESeq dataset
  }
}
