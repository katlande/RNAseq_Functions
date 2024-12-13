add_col <- function(df, colname, colval){
  df$v <- colval
  colnames(df)[ncol(df)] <- colname
  return(df)
}
