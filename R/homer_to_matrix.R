homer_to_matrix <- function (raw, cutchar = "_", instance = "1", pool = FALSE, round=TRUE) 
{
  raw <- raw[c(8:ncol(raw))]
  colnames(raw)[c(1)] <- c("Gene")
  if (instance == 1) {
    colnames(raw) <- gsub(paste0(cutchar, ".*"), "", colnames(raw))
  }
  else {
    colnames(raw) <- sub(sprintf(paste0("^((.*?", cutchar, 
                                        "){%d}.*?)", cutchar, ".*"), (instance - 1)), "\\1", 
                         colnames(raw))
  }
  raw$Gene <- sub("[|][|].*", "", raw$Gene)
  raw$Gene <- sub("[|].*", "", raw$Gene)
  if (pool == T) {
    raw <- aggregate(raw[c(2:ncol(raw))], by = list(raw$Gene), 
                     FUN = "sum")
    colnames(raw)[c(1)] <- c("Gene")
  }
  row.names(raw) <- raw$Gene
  raw <- raw[-c(1)]
  
  if(round==T){
    r <- as.data.frame(apply(raw, 2, as.integer))
    row.names(r) <- row.names(raw)
    raw <- r
  }
  
  return(raw)
}
