Reliable_heatmap <- function(norm, genes=NULL, sample_meta=NULL, annotation=F, annots=NULL, title="heatmap", 
         nrow=1, ncol=1, log=F, angle="45", frow=8, fcol=8, annots_as_factors=T, rowlabs=T, collabs=T,
         cols=c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D" , "#B2182B")){
  
  
  if(! is.null(genes)){
    n <- norm[row.names(norm)%in%genes,]
  }
 
  if(log==T){
    n <- log(n+1)
  }
  
  
  if(! is.null(sample_meta)){
    sample_meta <- sample_meta[row.names(sample_meta) %in% colnames(n),]
    n <- n[colnames(n) %in% row.names(sample_meta)]
    if(ncol(sample_meta) == 1){
       sample_meta$`Specify_annots_argument_to_remove_me` <- 1
    }  
  }
  
  if(annotation==T){
    
    if(is.null(sample_meta)){
      stop("Error: sample_meta must be provided if annotation is TRUE.")
    }
    
    if(is.null(annots)){
      s <- sample_meta
    } else {
      s <- sample_meta[colnames(sample_meta) %in% annots]
      if(annots_as_factors){
        s <- as.data.frame(apply(s, 2, as.factor))
      }
    }
    
    pheatmap::pheatmap(n, scale="row", color = cols, border_color = "black", 
                       clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
                       treeheight_row = 15, treeheight_col = 15, cutree_rows = nrow, cutree_cols = ncol,
                       fontsize_row = frow, fontsize_col = fcol, main = title, 
                       show_rownames = rowlabs, show_colnames = collabs,
                       annotation_col = s, angle_col = angle, annotation_names_col = F)->p
    
  } else{
    pheatmap::pheatmap(n, scale="row", color = cols, border_color = "black", 
                       clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
                       treeheight_row = 15, treeheight_col = 15, cutree_rows = nrow, cutree_cols = ncol,
                       fontsize_row = frow, fontsize_col = fcol,  main = title,
                       show_rownames = rowlabs, show_colnames = collabs) ->p
  }
  
  return(p)
}
