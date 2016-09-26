rm(list=ls())
X<-read.csv("./SampleData/FGD.csv",header=FALSE)
#For debugging of original LinkCluE
K<-4
trueLabels<-read.csv("./SampleData/FGT.csv",header=FALSE)
trueLabels<-trueLabels[,1]
M<-10
k<-ceiling(sqrt(nrow(X)))
scheme<-1
dcCTS<-0.8
dcSRS<-0.8
R<-5
dcASRS<-0.8
answer<-linkclueMain(X,M,k,scheme,K,dcCTS,dcSRS,R,dcASRS,trueLabels)
answer1<-linkclueMain(X,M,k,scheme,K,dcCTS,dcSRS,R,dcASRS,trueLabels)

# E<-crEnsemble(X,M,k,scheme)
# mNames<-c("cts","srs","asrs")
# dcs<-c(dcCTS,dcSRS,dcASRS)
# CR<-NULL
# for(i in 1:length(mNames)){
#   if(mNames[i]=="srs"){
#     S<-do.call(mNames[i],list(E,dcs[i],R))
#   }else{
#     S<-do.call(mNames[i],list(E,dcs[i]))
#   }
#   CR<-cbind(CR,clHC(S,K))
# }

#For use with ConClust and its outputs
library(magrittr)
library(tidyr)
library(reshape2)
library(biostatUtil)
library(tcgaHGSC)

data("hgsc")
dat <- hgsc[,-1]
rownames(dat) <- hgsc[,1]
tdat <- t(dat)
mdat <- melt(tdat) %>%
  set_colnames(c("patientID","Gene","Expr"))

E<-crEnsemble(X=tdat,M=NULL,k=4,scheme=3,reps=100)


#For debugging class relabelling
load("EforClassRelabel.rda")
k<-4
if(is.na(dim(E)[3])){
  E_imputed<-array(dim=c(dim(E),1))
}else{
  E_imputed<-array(dim=dim(E))
}
for(i in 1:dim(E_imputed)[3]){
  E_imputed[,,i]<-apply(E[,,i],2,knnImpute,tdat=X)
}
E_cs<-consensusSummaries(E_imputed,k=k)
for(i in 1:length(E_cs)){
  E_imputed[,,i]<-apply(E_imputed[,,i],2,classRelabel,cl.ref=E_cs[[i]]$consensusClass)
}
