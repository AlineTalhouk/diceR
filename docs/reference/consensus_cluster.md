# Consensus clustering

Runs consensus clustering across subsamples of the data, clustering
algorithms, and cluster sizes.

## Usage

``` r
consensus_cluster(
  data,
  nk = 2:4,
  p.item = 0.8,
  reps = 1000,
  algorithms = NULL,
  nmf.method = c("brunet", "lee"),
  hc.method = "average",
  xdim = NULL,
  ydim = NULL,
  rlen = 200,
  alpha = c(0.05, 0.01),
  minPts = 5,
  distance = "euclidean",
  abs = TRUE,
  prep.data = c("none", "full", "sampled"),
  scale = TRUE,
  type = c("conventional", "robust", "tsne"),
  min.var = 1,
  progress = TRUE,
  seed.nmf = 123456,
  seed.data = 1,
  file.name = NULL,
  time.saved = FALSE
)
```

## Arguments

- data:

  data matrix with rows as samples and columns as variables

- nk:

  number of clusters (k) requested; can specify a single integer or a
  range of integers to compute multiple k

- p.item:

  proportion of items to be used in subsampling within an algorithm

- reps:

  number of subsamples

- algorithms:

  vector of clustering algorithms for performing consensus clustering.
  Must be any number of the following: "nmf", "hc", "diana", "km",
  "pam", "ap", "sc", "gmm", "block", "som", "cmeans", "hdbscan". A
  custom clustering algorithm can be used.

- nmf.method:

  specify NMF-based algorithms to run. By default the "brunet" and "lee"
  algorithms are called. See
  [`NMF::nmf()`](https://rdrr.io/pkg/NMF/man/nmf.html) for details.

- hc.method:

  agglomeration method for hierarchical clustering. The the "average"
  method is used by default.
  See[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html) for
  details.

- xdim:

  x dimension of the SOM grid

- ydim:

  y dimension of the SOM grid

- rlen:

  the number of times the complete data set will be presented to the SOM
  network.

- alpha:

  SOM learning rate, a vector of two numbers indicating the amount of
  change. Default is to decline linearly from 0.05 to 0.01 over `rlen`
  updates. Not used for the batch algorithm.

- minPts:

  minimum size of clusters for HDBSCAN. Default is 5.

- distance:

  a vector of distance functions. Defaults to "euclidean". Other options
  are given in [`stats::dist()`](https://rdrr.io/r/stats/dist.html). A
  custom distance function can be used.

- abs:

  only used for `distance = c("spearman", "pearson")`. If `TRUE`, the
  absolute value is first applied to the distance before subtracting
  from 1, e.g., we use 1 - \|SCD\| instead of 1 - SCD for the spearman
  correlation distance.

- prep.data:

  Prepare the data on the "full" dataset, the "sampled" dataset, or
  "none" (default).

- scale:

  logical; should the data be centered and scaled?

- type:

  if we use "conventional" measures (default), then the mean and
  standard deviation are used for centering and scaling, respectively.
  If "robust" measures are specified, the median and median absolute
  deviation (MAD) are used. Alternatively, we can apply "tsne" for
  dimension reduction.

- min.var:

  minimum variability measure threshold used to filter the feature space
  for only highly variable features. Only features with a minimum
  variability measure across all samples greater than `min.var` will be
  used. If `type = "conventional"`, the standard deviation is the
  measure used, and if `type = "robust"`, the MAD is the measure used.

- progress:

  logical; should a progress bar be displayed?

- seed.nmf:

  random seed to use for NMF-based algorithms

- seed.data:

  seed to use to ensure each algorithm operates on the same set of
  subsamples

- file.name:

  if not `NULL`, the returned array will be saved at each iteration as
  well as at the end of the function call to an `rds` object with
  `file.name` as the file name.

- time.saved:

  logical; if `TRUE`, the date saved is appended to `file.name`. Only
  applicable when `file.name` is not `NULL`.

## Value

An array of dimension `nrow(x)` by `reps` by `length(algorithms)` by
`length(nk)`. Each cube of the array represents a different k. Each
slice of a cube is a matrix showing consensus clustering results for
algorithms. The matrices have a row for each sample, and a column for
each subsample. Each entry represents a class membership.

When "hdbscan" is part of `algorithms`, we do not include its clustering
array in the consensus result. Instead, we report two summary statistics
as attributes: the proportion of outliers and the number of clusters.

## Details

See examples for how to use custom algorithms and distance functions.
The default clustering algorithms provided are:

- "nmf": Nonnegative Matrix Factorization (using Kullback-Leibler
  Divergence or Euclidean distance; See Note for specifications.)

- "hc": Hierarchical Clustering

- "diana": DIvisive ANAlysis Clustering

- "km": K-Means Clustering

- "pam": Partition Around Medoids

- "ap": Affinity Propagation

- "sc": Spectral Clustering using Radial-Basis kernel function

- "gmm": Gaussian Mixture Model using Bayesian Information Criterion on
  EM algorithm

- "block": Biclustering using a latent block model

- "som": Self-Organizing Map (SOM) with Hierarchical Clustering

- "cmeans": Fuzzy C-Means Clustering

- "hdbscan": Hierarchical Density-based Spatial Clustering of
  Applications with Noise (HDBSCAN)

The progress bar increments on every unit of `reps`.

## Note

The `nmf.method` options are "brunet" (Kullback-Leibler Divergence) and
"lee" (Euclidean distance). When "hdbscan" is chosen as an algorithm to
use, its results are excluded from the rest of the consensus clusters.
This is because there is no guarantee that the cluster assignment will
have every sample clustered; more often than not there will be noise
points or outliers. In addition, the number of distinct clusters may not
even be equal to `nk`.

## Author

Derek Chiu, Aline Talhouk

## Examples

``` r
data(hgsc)
dat <- hgsc[1:100, 1:50]

# Custom distance function
manh <- function(x) {
  stats::dist(x, method = "manhattan")
}

# Custom clustering algorithm
agnes <- function(d, k) {
  return(as.integer(stats::cutree(cluster::agnes(d, diss = TRUE), k)))
}

assign("agnes", agnes, 1)

cc <- consensus_cluster(dat, reps = 6, algorithms = c("pam", "agnes"),
distance = c("euclidean", "manh"), progress = FALSE)
str(cc)
#>  int [1:100, 1:6, 1:4, 1:3] 1 1 NA NA NA 1 1 NA 2 NA ...
#>  - attr(*, "dimnames")=List of 4
#>   ..$ : chr [1:100] "TCGA.04.1331_PRO.C5" "TCGA.04.1332_MES.C1" "TCGA.04.1336_DIF.C4" "TCGA.04.1337_MES.C1" ...
#>   ..$ : chr [1:6] "R1" "R2" "R3" "R4" ...
#>   ..$ : chr [1:4] "PAM_Euclidean" "PAM_Manh" "AGNES_Euclidean" "AGNES_Manh"
#>   ..$ : chr [1:3] "2" "3" "4"
```
