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
E_imputed <- apply(E, 2:4, knn_impute, data = dat)
E_imputed2 <- imputeMissing(E, data = dat, imputeALL = TRUE)
dim(E_imputed2)

# Do one or all of the following

# Obtain a consensus:
## Compute co-association matrix and hierarchical cluster
CO <- consensus_summary(E_imputed)

CO_cl <- consensus_combine(E, element = "class")
CO_mat <- consensus_combine(E, element = "matrix")

## Compute Link-based similarity Matrix and cluster
S <- apply(E_imputed2, 3, LCE)
LCE_cl <- consensus_class(S[[1]]$CTS, 4)

## Reorder and k modes
kmodes_cl <- k_modes(E_imputed2)

## Reorder and majority vote
majvot_cl <- majority_voting(E_imputed2)

# Combine all results
Final <- cbind(CO_cl, LCE_cl, kmodes_cl, majvot_cl)

# Evaluate different clustering using internal, external (if applicable) and graphical method
consensus_evaluate(t(hgsc[-1, -1]), Final, CO_mat)   # TODO: modify `Final` and `CO_mat` so it resembles `E` data structure

# Identify the best number of clustering and test against the null hypothesis that the data was generated from a unimodal distribution
