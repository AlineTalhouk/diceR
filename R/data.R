#' Example data for imputation
#'
#' The data set is intended to be tested with \code{\link{knn_impute}}
#'
#' @format An array with 489 rows, 10 columns, 12 slices
"E_imputed"

#' High grade serous cancer data from TCGA used to classify subtypes
#'
#' There are 321 genes and 489 samples. The first column is a gene identifier.
#'
#' @format A data frame with 321 rows and 490 variables
"hgsc"

#' Output of diceR::ConClust using HGSC TCGA data with all 12 default algorithms and 50 repetitions
#'
#' The 1st dimension indicates number of samples (489), second dimension is number of repetitions (50), third dimension is number of algorithms
#'
#' @format An arrary of dimension 489 by 50 by 12 by 1
"hgsc"