
# Name Constants ----------------------------------------------------------

NALG <- "nmf"
DALG <- c("hc", "diana", "km", "pam")
OALG <- c("ap", "sc", "gmm", "block", "som", "cmeans", "hdbscan")
ALG_NAMES <- c(NALG, DALG, OALG)


# NMF-Based ---------------------------------------------------------------

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
  cl <- suppressWarnings(
    apcluster::apclusterK(
      apcluster::negDistMat, x, k, verbose = FALSE)@idx
  ) %>%
    dplyr::dense_rank() %>%
    purrr::set_names(rownames(x))
  if (length(cl)) cl else NA
}

#' Spectral Clustering (Radial-Basis Kernel)
#' @noRd
sc <- function(x, k) {
  kernlab::specc(as.matrix(x), k, kernel = "rbfdot")@.Data %>%
    purrr::set_names(rownames(x))
}

#' Gaussian Mixture Model
#' @noRd
gmm <- function(x, k) {
  mclust::Mclust(x, k, verbose = FALSE)$classification
}

#' Block Clustering (Co-clustering)
#' @noRd
block <- function(x, k) {
  cl <- tryCatch(
    blockcluster::cocluster(x, "continuous",
                            nbcocluster = c(k, k))@rowclass + 1,
    error = function(e) return(NA)
  )
  if (length(cl)) cl else NA
}

#' Self-Organizing Maps
#' @noRd
som <- function(x, k, xdim, ydim, rlen, alpha, method = "average") {
  if (!is.matrix(x)) x <- as.matrix(x)
  model <- som_train(x = x, xdim = xdim, ydim = ydim, rlen = rlen,
                     alpha = alpha)
  som_cluster(model = model, k = k, method = method)
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
#' @noRd
cmeans <- function(x, k) {
  fuzzy <- seq(1.1, 3, by = 0.1) %>%
    purrr::map(~ e1071::cmeans(x = x, centers = k, m = .x))
  mbest <- c("xie.beni", "fukuyama.sugeno", "partition.entropy") %>%
    purrr::map(
      ~ purrr::map(fuzzy,
                   function(f) e1071::fclustIndex(y = f, x = x, index = .x)
      )
    ) %>%
    purrr::map_int(which.min) %>%
    table() %>%
    which.max() %>%
    names() %>%
    as.integer()
  fuzzy[[mbest]]$cluster
}

#' Hierarchical Density-Based Spatial Clustering of Applications with Noise
#' summarize the proportion of outliers and number of clusters
#' remove from consensus array and assign as an attribute, if used
#' @noRd
hdbscan_summarize <- function(arr, algorithms) {
  if ("hdbscan" %in% algorithms) {
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
  }
  arr
}