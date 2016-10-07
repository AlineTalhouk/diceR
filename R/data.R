#' Example data for Linkage Cluster Ensemble (LCE)
#' 
#' The data set is intended to be used with the LCE family of functions
#' 
#' @format A matrix with 100 rows (samples) and 10 columns (cluster assignments).
"E_LCE"

#' Example data for imputation
#'
#' The data set is intended to be tested with \code{\link{knn_impute}}
#'
#' @format An array with 489 rows, 10 columns, 12 slices
"E_imputed"

#' Example dataset FGD
#' 
#' The data is simulated for cluster analysis.
#'
#' @format A data frame with 100 rows and 12 columns
"FGD"

#' Example dataset FGT
#'
#' The single column has 25 elements evenly placed in 4 groups.
#'
#' @format A data frame with 100 rows and 1 column
"FGT"

#' Example dataset LD
#'
#' The data is simulated for cluster analysis.
#'
#' @format A data frame with 72 rows and 1081 columns
"LD"

#' Example dataset LT
#'
#' The single column has 24 elements in group 1 and 48 elements in group 2.
#'
#' @format A data frame with 72 rows and 1 column
"LT"

#' Example dataset final_c1_valid_RandIndex
#'
#' This vector is intended to be tested in \code{ev_rand}
#'
#' @format A vector of length 200
"final_c1_valid_RandIndex"

#' Example dataset final_c2_valid_RandIndex
#'
#' This vector is intended to be tested in \code{ev_rand}
#'
#' @format A vector of length 200
"final_c2_valid_RandIndex"

#' TCGA example data
#'
#' The columns are intended to represent different cluster assignments.
#'
#' @format A matrix with 489 rows and 9 columns
"fixedKmeansTCGA_CR"

#' TCGA example data
#' 
#' The data is assessed using different similarity measures and external
#' validity measures.
#' 
#' @format A matrix with 3 rows and 9 columns
"fixedKmeansTCGA_V"

#' High grade serous cancer data from TCGA used to classify subtypes
#'
#' There are 321 genes and 489 samples. The first column is a gene identifier.
#'
#' @format A data frame with 321 rows and 490 variables
"hgsc"
