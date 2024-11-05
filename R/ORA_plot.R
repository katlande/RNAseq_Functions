ORA_plot <- function(ORA, balance=T, n=10, ontology_name="Biological Process"){
  
  if(balance == T){
    tmp <- rbind(ORA[ORA$enrichmentRatio > 0,][1:ceiling(n/2),],
                 ORA[ORA$enrichmentRatio < 0,][1:ceiling(n/2),])
  } else {
    tmp <- ORA[1:n,]
  }
  
  tmp$ID <- factor(tmp$ID, tmp$ID[order(tmp$enrichmentRatio)])
  tmp$plot_FDR <- -log10(tmp$FDR+min(tmp$FDR[!tmp$FDR==0]*0.5))
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
  
  ggplot(tmp, aes(y=ID, x=enrichmentRatio, 
                  size=plot_FDR, 
                  fill=sign(enrichmentRatio)*plot_FDR))+
    geom_point(aes(shape=sig))+
    scale_shape_manual("Significant", values = c(`FALSE`=13,`TRUE`=21), drop=F)+
    scale_fill_gradientn("pAdj", colours = rev(RColorBrewer::brewer.pal(9, "RdBu")), breaks=breaks, labels=labels, limits=c(-max_fdr, max_fdr))+
    scale_size_continuous(guide = "none")+
    xlab("Enrichment Ratio")+
    ylab("")+
    Ol_Reliable()+
    labs(title="Top ORA Hits", subtitle=ontology_name) -> plot
  return(plot)
  
}
