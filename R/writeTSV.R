writeTSV <- function(df, dir, rows=F){
  write.table(df, paste0(dir,".tsv"), quote=F, col.names=T, row.names=rows, sep="\t")
}
   
