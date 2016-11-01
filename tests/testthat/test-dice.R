library(dplyr)
data("hgsc")
dat <- t(hgsc[,-1])

exp.dat <- data.frame(initCol=rownames(dat)) %>%
    tidyr::separate(.,initCol,into=c("patientID","Class"), sep="_") 
refClass <- factor(exp.dat$Class)


c <- dice(dat, nk=4, algorithms = c("hcAEucl"), 
     consensusFUNS = c("kmodes", "majority", "LCE"))

c$clusters
