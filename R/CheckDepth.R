CheckDepth <- function(raw, minDepth=1e07){
  ggplot2::gplot(data.frame(nReads=colSums(raw), x="a"), ggplot2::aes(y=nReads, x=x))+
    ggplot2::geom_boxplot()+ ggplot2::geom_jitter(width = 0.2, height = 0)+
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank())+ ggplot2::ylab("Total Reads per Sample")+
    ggplot2::geom_hline(yintercept = minDepth, linetype="dashed", colour="red")+
    ggplot2::ggtitle("Read Depth")+ Ol_Reliable() -> g
  return(g)
}
