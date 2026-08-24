# Evaluate, trim, and reweigh algorithms

Evaluates algorithms on internal/external validation indices. Poor
performing algorithms can be trimmed from the ensemble. The remaining
algorithms can be given weights before use in consensus functions.

## Usage

``` r
consensus_evaluate(
  data,
  ...,
  cons.cl = NULL,
  ref.cl = NULL,
  k.method = NULL,
  plot = FALSE,
  trim = FALSE,
  reweigh = FALSE,
  n = 5,
  lower = 0,
  upper = 1
)
```

## Arguments

- data:

  data matrix with rows as samples and columns as variables

- ...:

  any number of objects outputted from
  [`consensus_cluster()`](https://alinetalhouk.github.io/diceR/reference/consensus_cluster.md)

- cons.cl:

  matrix of cluster assignments from consensus functions such as
  `kmodes` and `majority_voting`

- ref.cl:

  reference class

- k.method:

  determines the method to choose k when no reference class is given. If
  `ref.cl` is not `NULL`, this is the number of distinct classes in the
  reference; otherwise the chosen k is determined by the one giving the
  largest mean PAC across algorithms. Alternatively, specifying an
  integer will override the best chosen k, and specifying "all" will
  produce consensus results for all k values ("all" is implicitly used
  when there is only one k).

- plot:

  logical; if `TRUE`, `graph_all` is called

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

- lower:

  the lower bound that determines what is ambiguous

- upper:

  the upper bound that determines what is ambiguous

## Value

`consensus_evaluate` returns a list with the following elements

- `k`: if `ref.cl` is not `NULL`, this is the number of distinct classes
  in the reference; otherwise the chosen `k` is determined by the one
  giving the largest mean PAC across algorithms

- `pac`: a data frame showing the PAC for each combination of algorithm
  and cluster size

- `ii`: a list of data frames for all k showing internal evaluation
  indices

- `ei`: a data frame showing external evaluation indices for `k`

- `trim.obj`: A list with 4 elements

  - `alg.keep`: algorithms kept

  - `alg.remove`: algorithms removed

  - `rank.matrix`: a matrix of ranked algorithms for every internal
    evaluation index

  - `top.list`: final order of ranked algorithms

  - `E.new`: A new version of a `consensus_cluster` data object

## Details

This function always returns internal indices. If `ref.cl` is not
`NULL`, external indices are additionally shown. Relevant graphical
displays are also outputted. Algorithms are ranked across internal
indices using Rank Aggregation. Only the top `n` algorithms are kept,
the rest are trimmed.

## Examples

``` r
# Consensus clustering for multiple algorithms
set.seed(911)
x <- matrix(rnorm(500), ncol = 10)
CC <- consensus_cluster(x, nk = 3:4, reps = 10, algorithms = c("ap", "km"),
progress = FALSE)

# Evaluate algorithms on internal/external indices and trim algorithms:
# remove those ranking low on internal indices
set.seed(1)
ref.cl <- sample(1:4, 50, replace = TRUE)
z <- consensus_evaluate(x, CC, ref.cl = ref.cl, n = 1, trim = TRUE)
str(z, max.level = 2)
#> List of 5
#>  $ k       : int 4
#>  $ pac     :'data.frame':    2 obs. of  3 variables:
#>   ..$ k : chr [1:2] "3" "4"
#>   ..$ AP: num [1:2] 0.505 0.48
#>   ..$ KM: num [1:2] 0.514 0.498
#>  $ ii      :List of 2
#>   ..$ 3:'data.frame':    2 obs. of  16 variables:
#>   ..$ 4:'data.frame':    2 obs. of  16 variables:
#>  $ ei      :List of 1
#>   ..$ 4:'data.frame':    2 obs. of  19 variables:
#>  $ trim.obj:List of 5
#>   ..$ alg.keep   : chr "KM"
#>   ..$ alg.remove : chr "AP"
#>   ..$ rank.matrix:List of 1
#>   ..$ top.list   :List of 1
#>   ..$ E.new      :List of 1
```
