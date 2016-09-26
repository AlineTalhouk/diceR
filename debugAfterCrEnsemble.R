load("~/Documents/LinkCluE/debuggingDataBCCRC/X.rda")
X<-dat
rm(dat)
load("~/Documents/LinkCluE/debuggingDataBCCRC/E_conClust.rda")
R<-5
truelabels<-NULL
K<-4
mNames<-c("cts","srs","asrs")
dcs<-c(dcCTS=0.8,dcSRS=0.8,dcASRS=0.8)
algNames<-c("nmfDiv", "nmfEucl", "hcAEucl", "hcDianaEucl", "kmEucl",
            "kmSpear", "pamEucl", "pamSpear", "apEucl",
            "scRbf", "gmmBIC", "biclust")
methods<-c("CTS-SL","CTS-CL","CTS-AL","SRS-SL","SRS-CL","SRS-AL","ASRS-SL","ASRS-CL","ASRS-AL")
CR<-array(dim=c(nrow(X),length(methods),length(algNames)))
if(is.null(truelabels)){
  V<-array(dim=c(3,length(methods),length(algNames)))
}else{
  V<-array(dim=c(6,length(methods),length(algNames)))
}

for(j in 1:dim(E)[2]){
  tempCR<-NULL
  tempV<-NULL
  for(i in 1:length(mNames)){
    if(mNames[i]=="srs"){
      S<-do.call(mNames[i],list(breakE(E,algNames,j),dcs[i],R))
    }else{
      S<-do.call(mNames[i],list(breakE(E,algNames,j),dcs[i]))
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
  CR[,,j]<-tempCR
  V[,,j]<-tempV
}
