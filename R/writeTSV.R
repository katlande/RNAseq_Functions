writeTSV <- function(df, dir, rows=F){
  write.table(df, dir, quote=F, col.names=T, row.names==rows, sep="\t")
}
   
