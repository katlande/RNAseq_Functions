forceC2R <- function(df, col, keep=F, colname="KEEP"){
  if(keep ==T){
    df <- tibble::rownames_to_column(df, colname)
  }
  row.names(df) <- c()
  return(tibble::column_to_rownames(df, col))
}