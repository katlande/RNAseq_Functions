SetScores <- function (norm, goi, meta, out="data.frame") 
{
  
  if(!out%in%c("data.frame", "matrix")){
    warning(paste(out, "is not a recognized option for variable 'out'; must be one of 'data.frame' or 'matrix'"))
  } else {
    
    tdf <- t(apply(norm, 1, FUN = function(x) {
      as.numeric(scale(x))
    }))
    row.names(tdf) <- row.names(norm)
    colnames(tdf) <- colnames(norm)
    if (is.list(goi)) {
      for (i in 1:length(goi)) {
        if (i == 1) {
          
          if(out == "data.frame"){
            goidf <- data.frame(Sample = colnames(tdf), 
                                Score = rowMedians(t(tdf[row.names(tdf) %in% 
                                                           goi[[i]], ])), Set = names(goi)[i])
          } else{
            goidf <- data.frame(Sample = colnames(tdf), 
                                Set = rowMedians(t(tdf[row.names(tdf) %in% goi[[i]], ])))
            colnames(goidf)[2] <- names(goi)[i]
          }
        }
        else {
          if(out == "data.frame"){
            tmp <- data.frame(Sample = colnames(tdf), Score = rowMedians(t(tdf[row.names(tdf) %in% 
                                                                                 goi[[i]], ])), Set = names(goi)[i])
            goidf <- rbind(goidf, tmp)
          } else{
            tmp <- data.frame(Sample = colnames(tdf), 
                              Set = rowMedians(t(tdf[row.names(tdf) %in% goi[[i]], ])))
            colnames(tmp)[2] <- names(goi)[i]
            goidf <- merge(goidf, tmp, by="Sample")
          }
        }
      }
    }
    else {
      goidf <- data.frame(Sample = colnames(tdf), Score = rowMedians(t(tdf[row.names(tdf) %in% goi, ])))
    }
    score_df <- meta %>% rownames_to_column("Sample") %>% merge(., 
                                                                goidf, by = "Sample", all.x = F, all.y = T)
    return(score_df)
  }
}
