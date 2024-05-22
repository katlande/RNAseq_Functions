Check_Filter <- function(D, thresh=0){

  keep <- rowSums(counts(D)) >= thresh #T/F vect for row sum >= thresh
  message(table(keep)) #how many genes are filtered
  D <- D[keep,] #filter D

  temp_df <- data.frame(log10(rowSums(counts(D))+1)) #log transform total counts

  # Plotting and aesthetic functionality below:
  colnames(temp_df) <- "counts"
  plot <- ggplot2::ggplot(temp_df, aes(x=counts)) +
    ggplot2::geom_rect(mapping=NULL, xmin=0, xmax = log10(thresh+1), ymin=-1, ymax=1, fill="red", alpha=0.005)+
    ggplot2::geom_density() + ylab("Density\n")+
    ggplot2::scale_x_continuous(breaks=c(seq(0, max(temp_df$counts), length.out=8)),
                       labels=floor(10^(seq(0, max(temp_df$counts), length.out=8))),
                       limits=c(0, (max(temp_df$counts)*1.05)), expand=c(0,0))+
    xlab("\nRaw Counts")+ Ol_Reliable()

  if(thresh > 0){
    plot+ggplot2::labs(title=paste("Threshold =", thresh),
              subtitle=paste0(formatC(table(keep)[1]/sum(table(keep))*100, digits=3), "% of Genes Removed\n"))+
      ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5, size=14, face="bold"),
            plot.subtitle=ggplot2::element_text(hjust=0.5, size=9, face="italic")) -> plot
  }else{
    plot+ggplot2::labs(title="Threshold = 0",
              subtitle="0% of Genes Removed\n")+
      ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5, size=14, face="bold"),
            plot.subtitle=ggplot2::element_text(hjust=0.5, size=9, face="italic")) -> plot
  }
  return(plot)
}
