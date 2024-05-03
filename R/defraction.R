defraction <- function(raw){
  r <- row.names(raw)
  raw <- as.data.frame(apply(raw, 2, floor))
  row.names(raw) <- r
  return(raw)
}
