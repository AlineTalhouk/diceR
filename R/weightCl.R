#' Title function to compute weight for each pair of clusters using their shared members (Jaccard coefficient)
#'
#' @param E : N by M matrix of cluster ensemble
#'
#' @return wcl: an weighted cluster matrix
#' @export
#'
#' @examples
weightCl<-function(E){
  N=nrow(E)
  no_allcl<-max(max(E))
  pc<-matrix(rep(0,N*no_allcl),nrow=N)

  for(i in 1:N){
    pc[i,E[i,]]<-1
  }

  wcl<-matrix(rep(0,no_allcl^2),nrow=no_allcl)

  for(i in 1:(no_allcl-1)){
    for(j in (i+1):no_allcl){
      tmp<-sum(as.numeric((pc[,i]+pc[,j]))>0)
      if (tmp>0){
        wcl[i,j]<-sum(as.numeric((pc[,i]+pc[,j]))==2)/tmp
      }
    }
  }

  wcl<-wcl+t(wcl)
  return(wcl)
}

