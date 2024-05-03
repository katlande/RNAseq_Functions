SetScores <- function(norm, goi, meta){

  tdf <- t(apply(norm, 1, FUN=function(x){as.numeric(scale(x))}))
  row.names(tdf) <- row.names(norm)
  colnames(tdf) <- colnames(norm)

  if(is.list(goi)){
    for(i in 1:length(goi)){

      if(i == 1){
        goidf <- data.frame(Sample=colnames(tdf),
                            Score=rowMedians(t(tdf[row.names(tdf)%in%goi[[i]],])),
                            Set=names(goi)[i])
      } else{
        tmp <- data.frame(Sample=colnames(tdf),
                            Score=rowMedians(t(tdf[row.names(tdf)%in%goi[[i]],])),
                            Set=names(goi)[i])
        goidf <- rbind(goidf, tmp)
      }
    }
  } else{
    goidf <- data.frame(Sample=colnames(tdf), Score=rowMedians(t(tdf[row.names(tdf)%in%goi,])))
  }

  meta %>% rownames_to_column("Sample") %>%
    merge(. , goidf, by="Sample", all.x=F, all.y=T) -> score_df

  return(score_df)
}
