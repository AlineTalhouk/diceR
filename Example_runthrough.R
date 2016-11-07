# Load data
library(diceR)
data("hgsc")
head(hgsc)
dat <- t(hgsc[, -1])
k <- 4


# Fix the number of cluster or loop over several k options

# Generate Diverse Cluster Ensemble for a fixed k
E <- ConClust(dat, nc = k, reps = 10, 
              method = c("hcAEucl", "kmEucl", "scRbf", "apEucl", "gmmBIC"), 
              save = TRUE, file.name = "outputs/E")
dim(E)

# E <- load("outputs/E.rds")
# Impute Missing Values using KNN and majority vote

E_imputed <- impute_missing(E, dat)
E_knn <- E_imputed$knn
E_complete <- E_imputed$complete
dim(E_complete)

# Do one or all of the following

# Obtain a consensus:
## Compute co-association matrix and hierarchical cluster
CO <- consensus_summary(E_knn)

CO_cl <- consensus_combine(E, element = "class")
CO_mat <- consensus_combine(E, element = "matrix")

## Compute Link-based similarity Matrix and cluster
LCE_cl <- LCE(E_complete, data = dat, k)

## Reorder and k modes
kmodes_cl <- k_modes(E_complete)

## Reorder and majority vote
majvot_cl <- majority_voting(E_complete)

# Combine all results
cons_cl <- cbind(LCE = LCE_cl, KModes = kmodes_cl, MajVot = majvot_cl)

# Evaluate different clustering using internal, external (if applicable) and
# graphical method
consensus_evaluate(data = dat, k = 4, E, cons.cl = cons_cl, plot = FALSE)

# Identify the best number of clustering and test against the null hypothesis
# that the data was generated from a unimodal distribution

###################################################################################################################################################

#Test run of ConClust with tcga data on GSC with all default algorithms of ConClust and repetition 50
#GSC job number 1591697, completed


rm(list = ls())
library(NMF)
library(stats)
library(utils)
library(caret)
library(magrittr)
library(mclust)
library(kernlab)
library(apcluster)
library(doParallel)
library(parallel)
library(foreach)
library(readr)
library(dplyr)
library(plyr)
load('/home/gliu/documents/diceR/data/hgsc.rda')
source('/home/gliu/documents/diceR/R/utils.R')
source('/home/gliu/documents/diceR/R/ConClust.R')
dat <- hgsc[, -1]
rownames(dat) <- hgsc[, 1]
tdat <- t(dat)
result_conClustTest <- ConClust(x = tdat, nc = 4, reps = 50)
saveRDS(result_conClustTest, file = "/home/gliu/documents/conClustTest_out.rds")

################################################################################################################################################

#Test run of unsupervised evaluation
#GSC job 1607423, still running as of 10:45AM, Nov. 7, 2016

library(NMF)
library(clue)
library(stats)
library(utils)
library(caret)
library(magrittr)
library(mclust)
library(kernlab)
library(blockcluster)
library(apcluster)
library(doParallel)
library(parallel)
library(foreach)
library(klaR)
library(readr)
library(class)
load('/home/gliu/documents/diceR/data/hgsc.rda')
source('/home/gliu/documents/diceR/R/utils.R')
file.sources = list.files(
  path = '/home/gliu/documents/diceR/R',
  pattern = "*.R$",
  full.names = TRUE,
  ignore.case = TRUE
)
sapply(file.sources, source, .GlobalEnv)
data <- hgsc[, -1]
rownames(data) <- hgsc[, 1]
data <- t(data) #data  with rows are samples
rm(hgsc)
load('/home/gliu/documents/diceR/data/conClust_tcga.rda')


#Begin evaluation
consensusFUNS = c("kmodes", "majority", "cts", "asrs", "srs")
n <- dim(data)[1]
ncf <- length(consensusFUNS)
imp.obj <- impute_missing(conClust_tcga, data)
Ecomp <- imp.obj$complete
Final <- matrix(NA, nrow = n, ncol = ncf,
                dimnames = list(rownames(data), consensusFUNS))
for (i in 1:ncf) {
  cat(i)
  Final[, i] <- switch(consensusFUNS[i],
                       kmodes = k_modes(Ecomp),
                       majority = majority_voting(Ecomp),
                       CSPA = majority_voting(Ecomp),
                       cts = LCE(drop(Ecomp), k = 4, sim.mat = "cts"),
                       asrs = LCE(drop(Ecomp), k = 4, sim.mat = "asrs"),
                       srs = LCE(drop(Ecomp), k = 4, sim.mat = "srs"), 
  )
}

save(Final, file = "/home/gliu/documents/testEvaluation/testEvaluation_final.rda")
####################################################################################################################################################
