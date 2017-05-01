#' Gene expression data for High Grade Serous Carcinoma from TCGA
#'
#' There are 489 samples measured on 321 genes. Sample IDs are in the row names
#' and gene names are in the column names. This data set is used for clustering
#' HGSC into subtypes with prognostic significance. The cluster assignments
#' obtained by TCGA are indicated by the last six characters of each row name in
#' \code{hgsc}: \code{MES.C1}, \code{IMM.C2}, \code{DIF.C4}, and \code{PRO.C5}
#'
#' @format A data frame with 489 rows and 321 columns.
"hgsc"
