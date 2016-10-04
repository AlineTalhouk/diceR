#' Title function used for relabelling clusters in the ensemble 'E'
#'
#' @param E : N by M matrix of cluster ensemble
#'
#' @return newE: N by M matrix of relabelled ensemble, no_allcl: total number of clusters in the ensemble
#' @export
#'
#' @examples
relabelCl<-function(E){
  assertthat::assert_that(is.matrix(E))
  assertthat::assert_that(is.numeric(E))
  N<-nrow(E)
  M<-ncol(E)
  for(i in 1:N){
    for(j in 1:M){
      if(!checkPosInt(E[i,j])){
        stop("Error in relabelCl: one of the entries in the input matrix is not a positive integer.")
      }
    }
  }
  newE<-matrix(rep(0,N*M),nrow=N)
  ucl<-sort(unique(E[,1]))
  if(max(E[,1]!=length(ucl))==1){
    for (j in 1:length(ucl)){
      newE[which(E[,1]==ucl[j]),1]<-j
    }
  }
  for(i in 2:M){
    ucl<-sort(unique(E[,i]))
    prevCl<-length(sort(unique(c(newE[,1:(i-1)]))))
    for(j in 1:length(ucl)){
      newE[E[,i]==ucl[j],i]<-prevCl+j
    }
  }
  no_allcl<-max(max(newE))
  return(list(no_allcl=no_allcl,newE=newE))
}
