#' Title : function to compute validity scores for clustering results
#'
#' @param X : a dataset, rows correspond to observations, columns correspond to variables (excluding class labels!)
#' @param CR : a matrix of clustering resutls (excluding row and column headers !)
#' @param methods : a vector of legend strings for bar plots
#' @param truelabels : vector for true labels for each data point
#'
#' @return V: matrix of cluster validity scores
#' @export
#'
#' @examples
cleval<-function(X,CR,methods,truelabels){
  V<-NULL
  for (i in 1:ncol(CR)){
    cp<-valid_compactness(X,CR[,i])
    tempDB<-valid_DbDunn(X,CR[,i])
    db<-tempDB$DB
    dunn<-tempDB$Dunn
    if(nargs()==4){
      tempRI<-valid_RandIndex(CR[,i],truelabels)
      AR<-tempRI$AR
      RI<-tempRI$RI
      ca<-valid_CA(CR[,i],truelabels)
      V<-cbind(V,rbind(cp,db,dunn,AR,RI,ca))
    }else{
      V<-cbind(V,rbind(cp,db,dunn))
    }
  }
  if(nargs()==4){
    xlabel<-c("CP","DB","Dunn","AR","RI","CA")
  }else{
    xlabel<-c("CP","DB","Dunn")
  }
  #colnames(V)<-methods
  return(V)
}
