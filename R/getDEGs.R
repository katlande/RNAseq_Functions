getDEGs <- function(results, mode="count", alpha=0.05, lfc=0.5, direction="all"){
  
  if(mode=="count"){
    u <- nrow(results[results$log2FoldChange > abs(lfc) & results$padj <= alpha,])
    d <-  nrow(results[results$log2FoldChange < -1*abs(lfc) & results$padj <= alpha,])
    return(paste(u+d, "DEGs:", u, "upregulated;", d, "downregulated."))
  } else if(mode=="vector"){
    u <- results$Gene[results$log2FoldChange > abs(lfc) & results$padj <= alpha]
    d <-  results$Gene[results$log2FoldChange < -1*abs(lfc) & results$padj <= alpha]
    
    if(direction=="all"){
      return(c(u,d))
    } else if(direction=="up"){
      return(u)
    } else if(direction=="down"){
      return(d)
    } else{
      warning(paste0("direction '", direction, "' is not recognized. Use 'all', 'up', or 'down'."))
    }
    
  } else if(mode=="table"){
    u <- results[results$log2FoldChange > abs(lfc) & results$padj <= alpha,]
    d <-  results[results$log2FoldChange < -1*abs(lfc) & results$padj <= alpha,]
    
    if(direction=="all"){
      return(rbind(u,d))
    } else if(direction=="up"){
      return(u)
    } else if(direction=="down"){
      return(d)
    } else{
      warning(paste0("direction '", direction, "' is not recognized. Use 'all', 'up', or 'down'."))
    }
    
  } else{
    warning(paste0("Mode '", mode, "' is not recognized. Use 'count', 'table', or 'vector'."))
  }
  
}
