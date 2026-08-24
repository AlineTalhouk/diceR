# Diverse Clustering Ensemble

Runs consensus clustering across subsamples, algorithms, and number of
clusters (k).

## Usage

``` r
dice(
  data,
  nk,
  p.item = 0.8,
  reps = 10,
  algorithms = NULL,
  k.method = NULL,
  nmf.method = c("brunet", "lee"),
  hc.method = "average",
  distance = "euclidean",
  cons.funs = c("kmodes", "majority", "CSPA", "LCE", "LCA"),
  sim.mat = c("cts", "srs", "asrs"),
  prep.data = c("none", "full", "sampled"),
  min.var = 1,
  seed = 1,
  seed.data = 1,
  trim = FALSE,
  reweigh = FALSE,
  n = 5,
  evaluate = TRUE,
  plot = FALSE,
  ref.cl = NULL,
  progress = TRUE,
  verbose = TRUE
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

- k.method:

  determines the method to choose k when no reference class is given. If
  `ref.cl` is not `NULL`, this is the number of distinct classes in the
  reference; otherwise the chosen k is determined by the one giving the
  largest mean PAC across algorithms. Alternatively, specifying an
  integer will override the best chosen k, and specifying "all" will
  produce consensus results for all k values ("all" is implicitly used
  when there is only one k).

- nmf.method:

  specify NMF-based algorithms to run. By default the "brunet" and "lee"
  algorithms are called. See
  [`NMF::nmf()`](https://rdrr.io/pkg/NMF/man/nmf.html) for details.

- hc.method:

  agglomeration method for hierarchical clustering. The the "average"
  method is used by default.
  See[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html) for
  details.

- distance:

  a vector of distance functions. Defaults to "euclidean". Other options
  are given in [`stats::dist()`](https://rdrr.io/r/stats/dist.html). A
  custom distance function can be used.

- cons.funs:

  consensus functions to use. Current options are "kmodes" (k-modes),
  "majority" (majority voting), "CSPA" (Cluster-based Similarity
  Partitioning Algorithm), "LCE" (linkage clustering ensemble), "LCA"
  (latent class analysis)

- sim.mat:

  similarity matrix; choices are "cts", "srs", "asrs".

- prep.data:

  Prepare the data on the "full" dataset, the "sampled" dataset, or
  "none" (default).

- min.var:

  minimum variability measure threshold used to filter the feature space
  for only highly variable features. Only features with a minimum
  variability measure across all samples greater than `min.var` will be
  used. If `type = "conventional"`, the standard deviation is the
  measure used, and if `type = "robust"`, the MAD is the measure used.

- seed:

  random seed for knn imputation reproducibility

- seed.data:

  seed to use to ensure each algorithm operates on the same set of
  subsamples

- trim:

  logical; if `TRUE`, algorithms that score low on internal indices will
  be trimmed out

- reweigh:

  logical; if `TRUE`, after trimming out poor performing algorithms,
  each algorithm is reweighed depending on its internal indices.

- n:

  an integer specifying the top `n` algorithms to keep after trimming
  off the poor performing ones using Rank Aggregation. If the total
  number of algorithms is less than `n` no trimming is done.

- evaluate:

  logical; if `TRUE` (default), validity indices are returned. Internal
  validity indices are always computed. If `ref.cl` is not `NULL`, then
  external validity indices will also be computed.

- plot:

  logical; if `TRUE`, `graph_all` is called and a summary evaluation
  heatmap of ranked algorithms vs. internal validity indices is plotted
  as well.

- ref.cl:

  reference class

- progress:

  logical; should a progress bar be displayed?

- verbose:

  logical; if `TRUE`, console will print out messages for major tasks in
  the consensus clustering

## Value

A list with the following elements

- E:

  raw clustering ensemble object

- Eknn:

  clustering ensemble object with knn imputation used on `E`

- Ecomp:

  flattened ensemble object with remaining missing entries imputed by
  majority voting

- clusters:

  final clustering assignment from the diverse clustering ensemble
  method

- indices:

  if `evaluate = TRUE`, shows cluster evaluation indices; otherwise
  `NULL`

## Details

There are three ways to handle the input data before clustering via
argument `prep.data`. The default is to use the raw data as-is ("none").
Or, we can enact
[`prepare_data()`](https://alinetalhouk.github.io/diceR/reference/prepare_data.md)
on the full dataset ("full"), or the bootstrap sampled datasets
("sampled").

## Author

Aline Talhouk, Derek Chiu

## Examples

``` r
data(hgsc)
dat <- hgsc[1:100, 1:50]
ref.cl <- strsplit(rownames(dat), "_") |>
  purrr::map_chr(2) |>
  factor() |>
  as.integer()
dice.obj <- dice(dat, nk = 4, reps = 5, algorithms = "hc", cons.funs =
"kmodes", ref.cl = ref.cl, progress = FALSE, verbose = FALSE)
str(dice.obj, max.level = 2)
#> List of 5
#>  $ E       : int [1:100, 1:5, 1, 1] 1 1 NA NA NA 1 1 NA 1 NA ...
#>   ..- attr(*, "dimnames")=List of 4
#>  $ Eknn    : int [1:100, 1:5, 1, 1] 1 1 1 1 1 1 1 1 1 3 ...
#>   ..- attr(*, "dimnames")=List of 4
#>  $ Ecomp   : num [1:100, 1:5, 1] 1 1 1 1 1 1 1 1 1 3 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ clusters: int [1:100, 1:2] 4 3 1 3 3 4 1 3 2 4 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ indices :List of 5
#>   ..$ k   : int 4
#>   ..$ pac :'data.frame': 1 obs. of  2 variables:
#>   ..$ ii  :List of 1
#>   ..$ ei  :List of 1
#>   ..$ trim:List of 5
```
