data("hgsc")
dat <- t(hgsc[,-1])
dice(dat, nk=4, algorithms = c("hcAEucl", "kmEucl"), 
     consensusFUNS = c("kmodes", "majority", "LCE"))
