library(diceR)
data("hgsc")
head(hgsc)

# Fix the number of cluster or loop over several k options

# Generate Diverse Cluster Ensemble for a fixed k
E <-ConClust(t(hgsc[,-1]), nc = 2:4, reps = 100, 
             method = c("hcAEucl", "kmEucl", "scRbf","apEucl","gmmBIC"), 
             save = TRUE, file.name = "outputs/E")


# Impute Missing Values using KNN and majority vote
E_imputed <- apply(E,2:4, knn_impute, data= t(hgsc[,-1]))
E_imputed2 <- imputeMissing(E,data= t(hgsc[,-1]),imputeALL=TRUE)

# Do one or all of the following

# Obtain a consensus:
## Compute co-association matrix and hierarchical cluster
CO <- consensus_summary(E_imputed)

obj <- unlist(CO[[2]], recursive = FALSE)
comb.class <- as.data.frame(obj[grep("consensus_class",names(obj))])
cons.mat<- Reduce('+',obj[grep("consensus_matrix",names(obj))])/5

ccomp <- consensus_evaluate(t(hgsc[,-1]), cl.mat = data.matrix(tmp), cons.mat = cons.mat)

cweigh <- consensus_weigh(ccomp$internal, top = 12)
cm_algs <- consensus_matrix(ccomb_class)

cm_algs_w <- consensus_matrix(ccomb_class,
                              weights = cweigh$Weight[match(colnames(ccomb_class), 
                                                            cweigh$Algorithms)])
gplots::heatmap.2(cm_algs, dendrogram = "none", Rowv = TRUE, Colv = TRUE,
                  trace = "none", labRow = NA, labCol = NA)

str(cm_algs)

## Compute Link-based similarity Matrix and cluster
## Reorder and k modes
kmodes_cl <- apply(E_imputed2, 3, k_modes)

## Reorder and majority vote
majvot_cl <- apply(E_imputed2, 3, majority_voting)

# Evaluate different clustering using internal, external (if applicable) and graphical methods

# Identify the best number of clustering and test against the null hypothesis that the data was generated from a unimodal distribution
