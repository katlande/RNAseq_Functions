Heatmap_Grob <- function(mat, font_row=8, font_col=8, tilt=T, show_row=T, show_col=T, border="black",
                         clustDist_col="correlation", clustDist_row="correlation", clustMethod="ward.D2",
                         fill_colours=NULL, title=NULL, legend_title=NULL, scale_rows=T){
  
  if(scale_rows==T){
    tmp <- t(apply(mat, 1, FUN = function(x) {
      as.numeric(scale(x))
    }))
    row.names(tmp) <- row.names(mat)
    colnames(tmp) <- colnames(mat)
    tmp -> mat
  }
  mat <- as.data.frame(mat)
  
  tree_row <- cluster_mat(mat, distance = clustDist_row, method = clustMethod)
  tree_col <- cluster_mat(t(mat), distance = clustDist_col, method = clustMethod)
  col_order <- colnames(mat)[tree_col$order]
  row_order <- row.names(mat)[tree_row$order]
  
  gg <- reshape2::melt(data.matrix(mat))
  colnames(gg) <- c("Gene", "Sample", "Value")
  gg$Gene <- factor(as.character(gg$Gene), levels=row_order)
  gg$Sample <- factor(gg$Sample, levels=col_order)
  
  g <- ggplot2::ggplot(gg, aes(x=Sample, y=Gene, fill=Value))+
    ggplot2::geom_tile(colour=border)+
    Ol_Reliable()+ 
    ggplot2::ylab("")+ 
    ggplot2::xlab("")+
    ggplot2::scale_y_discrete(expand=c(0,0), position = "right")+
    ggplot2::scale_x_discrete(expand=c(0,0))+
    ggplot2::theme(plot.margin = unit(c(0.1,0.1,0.1,1), "cm"))
  
  
  if(is.null(legend_title)){
    legend_title <- ""
  } 
  
  if(is.null(fill_colours)){
    g <- g+ggplot2::scale_fill_gradient2(legend_title, high = "#B2182B", low="#2166AC", mid="#F7F7F7", midpoint = 0)
  } else{
    g <- g+ggplot2::scale_fill_gradientn(legend_title, colours = fill_colours)
  }
  
  if(show_row==T){
    g <- g+ggplot2::theme(axis.text.y = ggplot2::element_text(size=font_row, colour="black"))
  } else {
    g <- g+ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
  }
  
  if(show_col==T){
    
    if(tilt==T){
      g <- g+ggplot2::theme(axis.text.x = ggplot2::element_text(size=font_row, colour="black", angle=45, vjust=1, hjust=1))
    } else{
      g <- g+ggplot2::theme(axis.text.x = ggplot2::element_text(size=font_row, colour="black"))
    }
    
  } else {
    g <- g+ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
  }
  
  
  if(! is.null(title)){
    g <- g+ggtitle(title)
  }
  
  return(g)
}


