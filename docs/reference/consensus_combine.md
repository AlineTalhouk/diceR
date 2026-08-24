# Combine algorithms

Combines results for multiple objects from
[`consensus_cluster()`](https://alinetalhouk.github.io/diceR/reference/consensus_cluster.md)
and outputs either the consensus matrices or consensus classes for all
algorithms.

## Usage

``` r
consensus_combine(..., element = c("matrix", "class"))
```

## Arguments

- ...:

  any number of objects outputted from
  [`consensus_cluster()`](https://alinetalhouk.github.io/diceR/reference/consensus_cluster.md)

- element:

  either "matrix" or "class" to extract the consensus matrix or
  consensus class, respectively.

## Value

`consensus_combine` returns either a list of all consensus matrices or a
data frame showing all the consensus classes

## Details

This function is useful for collecting summaries because the original
results from `consensus_cluster` were combined to a single object. For
example, setting `element = "class"` returns a matrix of consensus
cluster assignments, which can be visualized as a consensus matrix
heatmap.

## Author

Derek Chiu

## Examples

``` r
# Consensus clustering for multiple algorithms
set.seed(911)
x <- matrix(rnorm(500), ncol = 10)
CC1 <- consensus_cluster(x, nk = 3:4, reps = 10, algorithms = "ap",
progress = FALSE)
CC2 <- consensus_cluster(x, nk = 3:4, reps = 10, algorithms = "km",
progress = FALSE)

# Combine and return either matrices or classes
y1 <- consensus_combine(CC1, CC2, element = "matrix")
str(y1)
#> List of 2
#>  $ 3:List of 2
#>   ..$ AP: num [1:50, 1:50] 1 0.429 0.167 0.429 1 ...
#>   ..$ KM: num [1:50, 1:50] 1 0.625 0 0.429 1 ...
#>  $ 4:List of 2
#>   ..$ AP: num [1:50, 1:50] 1 0.143 0.167 0 1 ...
#>   ..$ KM: num [1:50, 1:50] 1 0.25 0 0.125 1 ...
y2 <- consensus_combine(CC1, CC2, element = "class")
str(y2)
#> List of 2
#>  $ 3: int [1:50, 1:2] 1 1 2 1 1 1 2 1 3 2 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:2] "AP" "KM"
#>  $ 4: int [1:50, 1:2] 1 2 3 2 1 2 3 3 4 3 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:2] "AP" "KM"
```
