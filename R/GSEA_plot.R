GSEA_plot <- function(GSEA, balance=T, n=10, ontology_name="KEGG", logSize=F){
  
  if(balance == T){
    if(any(GSEA$FDR[GSEA$normalizedEnrichmentScore > 0][1:ceiling(n/2)] > 0.05) == T & any(GSEA$FDR[GSEA$normalizedEnrichmentScore < 0][1:ceiling(n/2)] > 0.05) == F){
  
  downlen <- length(which(GSEA$FDR[GSEA$normalizedEnrichmentScore < 0][1:ceiling(n)] < 0.05))
  uplen <- length(which(GSEA$FDR[GSEA$normalizedEnrichmentScore > 0][1:ceiling(n)] < 0.05))
  
  if(uplen+downlen > n){
    tmp <- rbind(GSEA[GSEA$normalizedEnrichmentScore > 0,][1:uplen,],
                 GSEA[GSEA$normalizedEnrichmentScore < 0,][1:(n-uplen),])
  } else {
    dif <- ceiling((n-(uplen+downlen))/2)
    tmp <- rbind(GSEA[GSEA$normalizedEnrichmentScore > 0,][1:(uplen+dif),],
                 GSEA[GSEA$normalizedEnrichmentScore < 0,][1:(downlen+dif),])
  }
  
  
} else if(any(GSEA$FDR[GSEA$normalizedEnrichmentScore > 0][1:ceiling(n/2)] > 0.05) == F & any(GSEA$FDR[GSEA$normalizedEnrichmentScore < 0][1:ceiling(n/2)] > 0.05) == T){
  downlen <- length(which(GSEA$FDR[GSEA$normalizedEnrichmentScore < 0][1:ceiling(n)] < 0.05))
  uplen <- length(which(GSEA$FDR[GSEA$normalizedEnrichmentScore > 0][1:ceiling(n)] < 0.05))
  
  if(uplen+downlen > n){
    tmp <- rbind(GSEA[GSEA$normalizedEnrichmentScore > 0,][1:(n-downlen),],
                 GSEA[GSEA$normalizedEnrichmentScore < 0,][1:downlen,])
  } else {
    dif <- ceiling((n-(uplen+downlen))/2)
    tmp <- rbind(GSEA[GSEA$normalizedEnrichmentScore > 0,][1:(uplen+dif),],
                 GSEA[GSEA$normalizedEnrichmentScore < 0,][1:(downlen+dif),])
  }
} else {
  tmp <- rbind(GSEA[GSEA$normalizedEnrichmentScore > 0,][1:ceiling(n/2),],
               GSEA[GSEA$normalizedEnrichmentScore < 0,][1:ceiling(n/2),])
}

  } else {
    tmp <- GSEA[1:n,]
  }
  
  tmp$ID <- factor(tmp$ID, tmp$ID[order(tmp$normalizedEnrichmentScore)])
  tmp$plot_FDR <- -log10(tmp$FDR+min(GSEA$FDR[!GSEA$FDR==0]*0.5))
  max_fdr <- ceiling(max(tmp$plot_FDR))
  tmp$sig <- factor(tmp$FDR <= 0.05, levels=c("TRUE", "FALSE"))
  
  if(max_fdr <= -log10(0.05)){
    max_fdr <- 2
    breaks <- c(log10(0.05), 0, -log10(0.05))
    labels <- c("0.05", "NS", "0.05")
  } else{
    breaks <- c(-1*floor(max(tmp$plot_FDR)), log10(0.05), -log10(0.05), floor(max(tmp$plot_FDR)))
    labels <- c((10^(-1*floor(max(tmp$plot_FDR)))), "0.05", "0.05",(10^(-1*floor(max(tmp$plot_FDR)))))
  }
  
  ggplot(tmp, aes(y=ID, x=normalizedEnrichmentScore, 
                  size=leadingEdgeNum, 
                  fill=sign(normalizedEnrichmentScore)*plot_FDR))+
    geom_point(aes(shape=sig))+
    #scale_shape_manual("Significant", values = c(`FALSE`=13,`TRUE`=21), drop=F)+
    scale_fill_gradientn("pAdj", colours = rev(RColorBrewer::brewer.pal(9, "RdBu")), breaks=breaks, labels=labels, limits=c(-max_fdr, max_fdr))+
    #scale_size_continuous(guide = "none")+
    xlab("Normalized Enrichment Score")+
    ylab("")+
    Ol_Reliable()+
    labs(title="Top GSEA Hits", subtitle=ontology_name) -> plot
  
  if(any(tmp$FDR > 0.05)){
    plot <- plot +  scale_shape_manual("Significant", values = c(`FALSE`=13,`TRUE`=21), drop=F)
  } else {
    plot <- plot +  scale_shape_manual("Significant", values = c(`FALSE`=13,`TRUE`=21), guide="none")
  }
  
  if(logSize){
    plot <- plot+ scale_size_continuous("nGenes", range = c(0.1, 5), trans = "log10")
  } else{
    plot <- plot+ scale_size_continuous("nGenes", range = c(0.1, 5))
  }
  
  
  return(plot)
}
