#' Title function to evaluate ensemble result using ConClust
#'
#' @param X : original data, rows are samples and columns are genes
#' @param E : output from crEnsemble in LinkCluE using scheme=3 (ConClust)
#' @param R : number of iterations for simrank algorithm (don't need for CTS)
#' @param truelabels :optional, known cluster labels for each data points
#' @param K :number of clusters in the consensus functions
#' @param algNames : names of algorithms used in ConClust
#' @param methods : default: "CTS-SL","CTS-CL","CTS-AL","SRS-SL","SRS-CL","SRS-AL","ASRS-SL","ASRS-CL","ASRS-AL"
#' @param mNames : methods for creating similarity matrices
#' @param dcs : decay constants vector
#'
#' @return
#' @export
#'
#' @examples
evalConCl<-function(X,E,R,truelabels=NULL,K,algNames=c("nmfDiv", "nmfEucl", "hcAEucl", "hcDianaEucl", "kmEucl",
                                                     "kmSpear", "pamEucl", "pamSpear", "apEucl",
                                                     "scRbf", "gmmBIC", "biclust"),
                    methods=c("CTS-SL","CTS-CL","CTS-AL","SRS-SL","SRS-CL","SRS-AL","ASRS-SL","ASRS-CL","ASRS-AL"),
                    mNames=c("cts","srs","asrs"),dcs=c(0.8,0.8,0.8)){
  print("Evaluating results from ConClust")
  CR<-array(dim=c(nrow(X),length(methods),length(algNames)))
  if(is.null(truelabels)){
    V<-array(dim=c(3,length(methods),length(algNames)))
    rownames(V)<-c("CP","DB","Dunn","AR","RI","CA")
  }else{
    V<-array(dim=c(6,length(methods),length(algNames)))
    rownames(V)<-c("CP","DB","Dunn")
  }
  for(l in 1:dim(E)[3]){
    tempCR<-NULL
    tempV<-NULL
    for(i in 1:length(mNames)){
      if(mNames[i]=="srs"){
        S<-do.call(mNames[i],list(data.matrix(E[,,l]),dcs[i],R))
      }else{
        S<-do.call(mNames[i],list(data.matrix(E[,,l]),dcs[i]))
      }
      tempCR<-cbind(tempCR,clHC(S,K))
    }
    if(is.null(truelabels)){
      tempV<-cleval(X,tempCR,methods)
    }else{
      tempV<-cleval(X,tempCR,methods,truelabels)
    }
    colnames(tempV)<-methods
    colnames(tempCR)<-methods
    CR[,,l]<-tempCR
    V[,,l]<-tempV
  }
  colnames(CR)<-methods
  colnames(V)<-methods
  rownames(CR)<-rownames(X)
  return(list(CR=CR,V=V))
}
