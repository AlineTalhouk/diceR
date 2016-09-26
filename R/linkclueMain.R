#' Title function to perform link-based cluster ensemble for a dataset
#'
#' @param X : a dataset whose rows are observations and columns are variables
#' @param M : preferred number of base clusterings in the ensemble
#' @param k : preferred number of clusters in the base clusterings
#' @param scheme : cluster gensemble generating scheme (1=fixed, 2=random k)
#' @param K : number of clusters in the consensus functions
#' @param dcCTS : decay factor for CTS method
#' @param dcSRS ; decay factor for SRS method
#' @param R : number of iterations for simrank algorithm (don't need for CTS)
#' @param dcASRS : decay factor for ASRS method
#' @param truelabels : optional, known cluster labels for each data points
#' @param reps: number of repetitions if scheme is 3 for conclust clustering
#'
#' @return
#' @export
#'
#' @examples
linkclueMain<-function(X,M=NULL,k,scheme,K,dcCTS=0.8,dcSRS=0.8,R=5,dcASRS=0.8,truelabels=NULL,reps=1000){
  mNames<-c("cts","srs","asrs")
  dcs<-c(dcCTS,dcSRS,dcASRS)
  methods<-c("CTS-SL","CTS-CL","CTS-AL","SRS-SL","SRS-CL","SRS-AL","ASRS-SL","ASRS-CL","ASRS-AL")
  algNames<-c("nmfDiv", "nmfEucl", "hcAEucl", "hcDianaEucl", "kmEucl",
                        "kmSpear", "pamEucl", "pamSpear", "apEucl",
                        "scRbf", "gmmBIC", "biclust")
  if(!any(c(1,2,3)==scheme)){
    stop("Scheme should be 1 or 2 for fixed k or random k or 3 for conclust")
  }
  if(any(dcs>1)|any(dcs<0)){
    stop("decay constants must be between 0 and 1")
  }
  # intVals<-c(M,k,K,R)
  # for(i in 1:length(intVals)){
  #   if(!checkPosInt(intVals[i])){
  #     stop("At least one of M,k,K, or R is not an integer")
  #   }
  # }
  #set.seed(1979)
  E<-crEnsemble(X,M,k,scheme,reps=reps) # n x M matrix of clusterings generated in crEnsemble
  if(scheme!=3){
    CR<-NULL
    for(i in 1:length(mNames)){
      if(mNames[i]=="srs"){
        S<-do.call(mNames[i],list(E,dcs[i],R))
      }else{
        S<-do.call(mNames[i],list(E,dcs[i]))
      }
      CR<-cbind(CR,clHC(S,K))
    }
    if(is.null(truelabels)){
      V<-cleval(X,CR,methods)
    }else{
      V<-cleval(X,CR,methods,truelabels)
    }
    colnames(V)<-methods
    colnames(CR)<-methods
    barplot(t(V),beside=TRUE,col=c("red","green","blue","grey","yellow","orange","purple","brown","black"),legend=rownames(t(V)),main = "Cluster Validity Comparison")
    return(list(CR=CR,V=V))
  }else{
    return(evalConCl(X=X,E=E,R=R,truelabels=truelabels,K=K))
  }
}
