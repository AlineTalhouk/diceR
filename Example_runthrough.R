library(diceR)
data("hgsc")
head(hgsc)

# Fix the number of cluster or loop over several k options
k=4
B=10

# Generate Diverse Cluster Ensemble for a fixed k
E <- ConClust(t(hgsc[,-1]), k = k, reps = B, 
              method = c("hcAEucl","kmEucl","scRbf"))

# Impute Missing Values using KNN and majority vote
E_imputed <- imputeMissing(E,data= t(hgsc[,-1]),imputeALL=TRUE)

# Do one or all of the following

# Obtain a consensus:
## Compute co-association matrix and hierarchical cluster
CO <- consensus_summary(E_imputed,k)
CO_cl <-consensus_class(consensus_combine(CO, element="class"), k = 4)

## Compute Link-based similarity Matrix and cluster
LCE <- link_clust(E = E_imputed, 
                     dcCTS = 0.8, dcSRS = 0.8, dcASRS = 0.8, R = 10,
                  is.relabelled = TRUE)
CTS_cl <- consensus_class(LCE)
## Reorder and k modes
kmodes_cl <- k_modes(E_imputed)

## Reorder and majority vote
majvot_cl <- majority_voting(E_imputed)

# Evaluate different clustering using internal, external (if applicable) and graphical methods
CCO <- consensus_combine(data.frame(CO_cl,kmodes_cl,majvot_cl), element = "class")
consensus_evaluate(data= t(hgsc[,-1]),data.frame(CO_cl,kmodes_cl,majvot_cl),cons.mat = CO)

# Identify the best number of clustering and test against the null hypothesis that the data was generated from a unimodal distribution
