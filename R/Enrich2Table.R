Enrich2Table <- function(df, mode="ORA", n=10){
  
  if(mode == "ORA"){
    t <- setNames(ORA[c(1,8,10,6,5)], c("Term", "Enrichment", "FDR", "nGenes", "Total Genes"))
  } else{
    t <- setNames(GSEA[c(1,6,8,11,9)], c("Term", "Enrichment", "FDR", "nGenes", "Total Genes"))
  }
  t$Enrichment <- formatC(t$Enrichment, digits = 2)
  t$FDR <- formatC(t$FDR, digits = 2, format = "E")
  
  g <- gridExtra::tableGrob(t[1:n,], rows = NULL)
  g <- gtable::gtable_add_grob(g, grobs = grid::rectGrob(gp = grid::gpar(fill = NA, lwd = 2)),
                               t = 2, b = nrow(g), l = 1, r = ncol(g))
  g <- gtable::gtable_add_grob(g, grobs = grid::rectGrob(gp = grid::gpar(fill = NA, lwd = 2)),
                               t = 1, l = 1, r = ncol(g))
  return(grid::grid.draw(g))
}
