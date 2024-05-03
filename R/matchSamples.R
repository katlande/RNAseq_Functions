matchSamples <- function(raw, filt.meta){
  filt.raw <- raw[colnames(raw) %in% row.names(filt.meta)]
  return(filt.raw)
} # match raw to meta
