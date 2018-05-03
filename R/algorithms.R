
# Name Constants ----------------------------------------------------------

NALG <- "nmf"
DALG <- c("hc", "diana", "km", "pam")
OALG <- c("ap", "sc", "gmm", "block", "som", "cmeans", "hdbscan")
ALG_NAMES <- c(NALG, DALG, OALG)


# NMF-Based ---------------------------------------------------------------

#' Nonnegative Matrix Factorization
#' Transpose since input for NMF::nmf uses rows as vars, cols as samples
#' @noRd
nmf <- function(x, k, method, seed) {
  NMF::predict(NMF::nmf(t(x), rank = k, method = method, seed = seed))
}

#' Transform to non-negative matrix by column-binding a negative replicate and
#' then coercing all negative values to 0
#' @noRd
nmf_transform <- function(x) {
  x %>%
    cbind(-.) %>%
    as.data.frame() %>%
    purrr::map_dfc(~ ifelse(.x < 0, 0, .x))
}

# Distance-Based ----------------------------------------------------------

#' Hierarchical Clustering
#' @noRd
hc <- function(d, k, method = "average") {
  as.integer(stats::cutree(stats::hclust(d, method = method), k))
}

#' DIvisive ANAlysis Clustering
#' @noRd
diana <- function(d, k) {
  as.integer(stats::cutree(cluster::diana(d, diss = TRUE), k))
}

#' K-Means Clustering
#' @noRd
km <- function(d, k) {
  as.integer(stats::kmeans(d, k)$cluster)
}

#' Partitioning Around Medoids
#' @noRd
pam <- function(d, k) {
  as.integer(cluster::pam(d, k, cluster.only = TRUE))
}


# All Other Algorithms ----------------------------------------------------

#' Affinity Propagation
#' @noRd
ap <- function(x, k) {
  suppressWarnings(
    apcluster::apclusterK(apcluster::negDistMat, x, k, verbose = FALSE)@idx
  ) %>%
    dplyr::dense_rank() %>%
    purrr::when(length(.) > 0 ~ ., ~ NA)
}

#' Spectral Clustering (Radial-Basis Kernel)
#' @noRd
sc <- function(x, k) {
  kernlab::specc(as.matrix(x), k, kernel = "rbfdot")@.Data
}

#' Gaussian Mixture Model
#' @noRd
gmm <- function(x, k) {
  mclust::Mclust(x, k, verbose = FALSE)$classification
}

#' Block Clustering (Co-clustering)
#' @noRd
block <- function(x, k) {
  tryCatch(
    sink_output(
      blockcluster::cocluster(as.matrix(x), "continuous",
                              nbcocluster = c(k, k))@rowclass + 1
    ),
    error = function(e) return(NA)
  ) %>%
    purrr::when(length(.) > 0 ~ ., ~ NA)
}

#' Self-Organizing Maps
#' @noRd
som <- function(x, k, xdim, ydim, rlen, alpha, method = "average") {
  x %>%
    purrr::when(is.matrix(.) ~ ., ~ as.matrix(.)) %>%
    som_train(xdim = xdim, ydim = ydim, rlen = rlen, alpha = alpha) %>%
    som_cluster(k = k, method = method)
}

#' Train the SOM, specifiy grid size and other optional parameters based on the
#' SOM Toolbox
#' @references
#'   http://www.cis.hut.fi/somtoolbox/package/docs2/som_topol_struct.html
#' @noRd
som_train <- function(x, xdim, ydim, rlen, alpha, topo = "hexagonal") {
  # Create SOM grid and map data into the grid
  neurons <- 5 * sqrt(nrow(x))
  eigenvalues <- eigen(stats::cor(x))$values
  eigenratio <- eigenvalues[1] / eigenvalues[2]
  xdim <- xdim %||% sqrt(neurons / eigenratio)
  ydim <- ydim %||% neurons / xdim
  grid <- kohonen::somgrid(xdim = xdim, ydim = ydim, topo = topo)
  kohonen::som(x, grid = grid, rlen = rlen, alpha = alpha)
}

#' Perform hc: cut tree into k groups and output cluster labels for original data
#' @noRd
som_cluster <- function(model, k, method) {
  # Get distance matrix, use hc to cluster the codebook vectors
  cl <- hc(stats::dist(kohonen::getCodes(model, 1)),
           k = k,
           method = method)
  pred <- stats::predict(model)$unit.classif
  cl[pred]
}

#' Fuzzy C-Means (using best m via validity/performance measures)
#'
#' Fuzzifier m is a function of N and D (Equation 5 from ref.). If centers
#' can't be calculated, use frequently used value of m = 2.
#'
#' @references https://academic.oup.com/bioinformatics/article/26/22/2841/227572
#' @noRd
cmeans <- function(x, k) {
  N <- nrow(x)
  D <- ncol(x)
  m <- 1 + (1418 / N + 22.05) * D ^ (-2) +
    (12.33 / N + 0.243) * D ^ (-0.0406 * log(N) - 0.1134)
  fuzzy <- e1071::cmeans(x = x, centers = k, m = m)
  if (length(fuzzy$cluster) == 0) {
    fuzzy <- e1071::cmeans(x = x, centers = k, m = 2)
  }
  fuzzy$cluster
}

#' Hierarchical Density-Based Spatial Clustering of Applications with Noise
#' @noRd
hdbscan <- function(x, minPts) {
  dbscan::hdbscan(x = x, minPts = minPts)$cluster
}

#' Summarize the proportion of outliers and number of clusters
#' remove from consensus array and assign as an attribute, if used
#' @noRd
hdbscan_summarize <- function(arr) {
  h.idx <- match("HDBSCAN", dimnames(arr)[[3]])
  h.obj <- arr[, , h.idx, ] %>%
    as.data.frame() %>%
    purrr::map(~ {
      c(prop_outlier = sum(.x == 0, na.rm = TRUE) / sum(!is.na(.x)),
        num_cluster = dplyr::n_distinct(!.x %in% c(NA, 0)))
    }) %>%
    purrr::transpose() %>%
    purrr::map(unlist)
  arr <- arr[, , -h.idx, , drop = FALSE]
  attr(arr, "hdbscan") <- h.obj
  arr
}
