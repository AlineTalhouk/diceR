data("hgsc")
dat <- hgsc[,-1]
rownames(dat) <- hgsc[,1]
tdat <- t(dat)
mdat <- melt(tdat) %>%
  set_colnames(c("patientID","Gene","Expr"))

E<-ClusterDuck::crEnsemble(X=tdat,M=NULL,k=4,reps=10)
E_compressed<-E[,,1]
if(!is.na(dim(E)[3])){
  for(i in 2:dim(E)[3]){
    E_compressed<-cbind(E_compressed,E[,,i])
  }
}
