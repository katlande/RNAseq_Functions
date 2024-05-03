makeComp <- function(Dnorm, var, up, down){
  df <- na.omit(as.data.frame(DESeq2::results(Dnorm, contrast = c(var, up, down))))
  df$Gene <- row.names(df)
  df$Dir <- down
  df$Dir[df$log2FoldChange > 0] <- up
  return(df)
}
