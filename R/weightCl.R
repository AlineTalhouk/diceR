#' Title function to compute weight for each pair of clusters using their shared members (Jaccard coefficient)
#'
#' @param E : N by M matrix of cluster ensemble
#'
#' @return wcl: an p by p weighted cluster matrix where p denotes number of classes
#' @export
#'
#' @examples
#' data("E_LCE")
#' wl<-weightCl(E_LCE)
weightCl<-function(E){
  assertthat::assert_that(is.matrix(E))
  assertthat::assert_that(is.numeric(E))
  for(i in 1:nrow(E)){
    for(j in 1:ncol(E)){
      if(!checkPosInt(E[i,j])){
        stop("Error in relabelCl: one of the entries in the input matrix is not a positive integer.")
      }
    }
  }
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

