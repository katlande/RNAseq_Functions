SetEnrichment <- function(geneList, bkgd_subset, bkgd){
  length(geneList[geneList%in%bkgd_subset]) ->a
  length(bkgd_subset[!bkgd_subset%in%geneList]) ->b
  length(geneList[!geneList%in%bkgd_subset])->c
  length(bkgd[!bkgd%in%c(bkgd_subset,geneList)])->d

  fisher.test(matrix(c(a,c,b,d), nrow = 2, ncol=2))$p.value->p
  (a/c)/(b/d)->e

  if(is.nan(e)){
    if((c==0 | d==0) & ! (a==0 | b==0)){
      e <- Inf
    } else{
      e<-0
    }
  }

  #return(c(p,e))
  return(c(p,e))
}
