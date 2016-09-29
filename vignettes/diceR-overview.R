## ---- echo=FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE,
                      comment = "#>")

## ------------------------------------------------------------------------
# install.packages("devtools")
# devtools::install_github("AlineTalhouk/diceR")
library(diceR)
data(hgsc)

## ---- results='hide'-----------------------------------------------------
dat <- hgsc[, -1]
CC <- ConClust(dat, k = 4, reps = 10, pItem = 0.8,
               method = c("hcAEucl", "kmSpear", "pamSpear", "biclust"), save = FALSE)

## ------------------------------------------------------------------------
str(CC)
knitr::kable(head(CC[, , "biclust"]))

## ------------------------------------------------------------------------
CC.impute <- apply(CC, 2:3, knn_impute, t(dat))
sum(is.na(CC))
sum(is.na(CC.impute))

## ------------------------------------------------------------------------
library(gplots)

# Consensus matrix for PAM only
cm <- consensus_matrix(CC[, , "pamSpear"])
heatmap.2(cm, dendrogram = "none", Rowv = TRUE, Colv = TRUE,
          trace = "none", labRow = NA, labCol = NA, key = FALSE)

