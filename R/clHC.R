#' Title: function to perform final clustering using hierarchical algorithms (single linkage-SL, complete linkage-CL, average linkage-AL)
#'
#' @param S : N by N similarity matrix
#' @param K : preferred number of clusters
#'
#' @return CR: an N by 3 matrix of clustering results from SL, CL, and AL
#' @export
#'
#' @examples
clHC<-function(S,K){
  CR<-NULL
  d<-LinkCluE::stod(S)
  Z<-LinkCluE::linkage(d,"single")
  CR<-cbind(CR,LinkCluE::cluster(Z,K))
  Z<-LinkCluE::linkage(d,"complete")
  CR<-cbind(CR,LinkCluE::cluster(Z,K))
  Z<-LinkCluE::linkage(d,"average")
  CR<-cbind(CR,LinkCluE::cluster(Z,K))
  return(CR)
}
