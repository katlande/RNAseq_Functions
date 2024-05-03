customPCA <- function(D, meta, ptcol, filter_threshold=5,  ptsize=2,
                      labsize=3, ptshape=NA, main="PCA", label=F, facet=NA,
                      frows=NA, fcols=NA, fscale="fixed", custom_cols=NULL){

  norm <- filt.norm(D, filter_threshold)

  if(! identical(row.names(meta), colnames(norm))){
    meta <- meta[row.names(meta) %in% colnames(norm),]
    norm <- norm[colnames(norm) %in% row.names(meta)]
  }

  pca <- prcomp(t(norm))
  varPC1 <- paste0(substr(as.character(((pca$sdev^2/sum(pca$sdev^2))*100)[1]),1,5),"%")
  varPC2 <- paste0(substr(as.character(((pca$sdev^2/sum(pca$sdev^2))*100)[2]),1,5),"%")
  temp_df <- as.data.frame(pca$x[])
  temp_df$sample <- row.names(temp_df)

  if(is.na(ptshape)){
    plot <- ggplot2::ggplot(data = temp_df, ggplot2::aes(x=PC1, y=PC2, color=as.factor(meta[[ptcol]])))+
      ggplot2::geom_jitter(size=ptsize, alpha=0.9)+
      ggplot2::labs(colour=ptcol)
  } else {
    plot <- ggplot2::ggplot(data = temp_df, ggplot2::aes(x=PC1, y=PC2, color=as.factor(meta[[ptcol]]),
                                                         shape=as.factor(meta[[ptshape]])))+
      ggplot2::geom_jitter(size=ptsize, alpha=0.9)+
      ggplot2::labs(colour=ptcol, shape=ptshape)
  }

  if(label==T){
    plot <- plot+ggrepel::geom_text_repel(size=labsize, label=temp_df$sample)
  }

  plot <- plot+
    ggplot2::ylab(paste0("PC2 - ", varPC2," of Variance"))+
    ggplot2::xlab(paste0("PC1 - ", varPC1," of Variance"))+
    ggplot2::ggtitle(main)+
    Ol_Reliable()


  if(! is.na(facet)){
    if(is.na(frows) & is.na(fcols)){
      plot <- plot+ggplot2::facet_wrap(~as.factor(meta[[facet]]), scales = fscale)
    } else if(! is.na(frows) & is.na(fcols)){
      plot <- plot+ggplot2::facet_wrap(~as.factor(meta[[facet]]), nrow=frows, scales = fscale)
    } else if(! is.na(fcols) & is.na(frows)){
      plot <- plot+ggplot2::facet_wrap(~as.factor(meta[[facet]]), ncol=fcols, scales = fscale)
    } else{
      plot <- plot+ggplot2::facet_wrap(~as.factor(meta[[facet]]), ncol=fcols, nrow=frows, scales = fscale)
    }
  }


  if(! is.null(custom_cols)){
    plot <- plot+ggplot2::scale_colour_manual(values=custom_cols)
  }

  return(plot)
}
