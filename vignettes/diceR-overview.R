## ----setup, echo=FALSE---------------------------------------------------
knitr::opts_chunk$set(
	message = FALSE,
	warning = FALSE,
	collapse = TRUE,
	comment = "#>",
	fig.align = "center",
	fig.width = 6,
	fig.height = 4.5
)

## ----load----------------------------------------------------------------
# install.packages("devtools")
# devtools::install_github("AlineTalhouk/diceR")
library(diceR)
library(dplyr)
library(ggplot2)
library(knitr)
data(hgsc)

## ----ConClust, results='hide'--------------------------------------------
dat <- t(hgsc[, -1])
CC <- ConClust(dat, k = 4, reps = 10, pItem = 0.8,
               method = c("hcAEucl", "kmSpear", "pamSpear", "biclust"), save = FALSE)

## ----ConClust_biclust----------------------------------------------------
str(CC)
kable(head(CC[, , "biclust"]))

## ----knn_impute----------------------------------------------------------
CC.impute <- apply(CC, 2:3, knn_impute, dat)
sum(is.na(CC))
sum(is.na(CC.impute))

## ----consensus_matrix----------------------------------------------------
cm <- consensus_matrix(CC[, , "pamSpear"])
gplots::heatmap.2(cm, dendrogram = "none", Rowv = TRUE, Colv = TRUE,
                  trace = "none", labRow = NA, labCol = NA)

## ----consensus_class-----------------------------------------------------
cclass <- consensus_class(cm, k = 4)
kable(as.data.frame(table(cclass)))

## ----consensus_summary, results='hide'-----------------------------------
cs <- consensus_summary(CC, k = 4)

## ----consensus_summary_str-----------------------------------------------
str(cs, max.level = 2)

## ----consensus_combine---------------------------------------------------
ccomb_matrix <- consensus_combine(cs, element = "matrix")
ccomb_class <- consensus_combine(cs, element = "class")

## ----consensus_combine_str-----------------------------------------------
str(ccomb_matrix, max.level = 2)
kable(head(ccomb_class))

## ----consensus_combine_2, results='hide'---------------------------------
CC2 <- ConClust(dat, k = 4, reps = 10, pItem = 0.8,
                method = c("hcDianaEucl", "pamEucl"), save = FALSE)
cs2 <- consensus_summary(CC2, k = 4)
ccomb_class2 <- consensus_combine(cs, cs2, element = "class")

## ----consensus_combine_2_str---------------------------------------------
kable(head(ccomb_class2))

## ----consensus_evaluate--------------------------------------------------
ccomp <- consensus_evaluate(dat, cl.mat = ccomb_class, cons.mat = ccomb_matrix)
kable(head(ccomp))

## ----consensus_weigh-----------------------------------------------------
cweigh <- consensus_weigh(ccomp$internal)
kable(head(cweigh))

## ----consensus_weigh_example---------------------------------------------
# Generate both consensus matrices
cm_algs <- consensus_matrix(ccomb_class)
cm_algs_w <- consensus_matrix(ccomb_class,
                              weights = cweigh$Weight[match(colnames(ccomb_class), cweigh$Algorithms)])

# Calculate range of PAC from 0 to Upper Limit
pac_dat <- data.frame(Upper = seq(0, 1, 0.01)) %>% 
  mutate(Unweighted = sapply(Upper, PAC, cm = cm_algs, lower = 0),
         Weighted = sapply(Upper, PAC, cm = cm_algs_w, lower = 0)) %>% 
  tidyr::gather(key = Type, value = PAC, 2:3)

# Show distribution of PAC for both Types of consensus matrices
ggplot(pac_dat, aes(Upper, PAC)) + 
  geom_line(aes(colour = Type)) + 
  labs(x = "Upper Limit",
       title = "Distribution of PAC from 0 to Upper Limit")

## ----consensus_confmat---------------------------------------------------
set.seed(10)
cclass.true <- sample(1:4, size = length(cclass), replace = TRUE)
table(cclass, cclass.true)
consensus_confmat(cclass, cclass.true)

