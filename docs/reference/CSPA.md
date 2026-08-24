# Cluster-based Similarity Partitioning Algorithm (CSPA)

Performs hierarchical clustering on a stack of consensus matrices to
obtain consensus class labels.

## Usage

``` r
CSPA(E, k)
```

## Arguments

- E:

  is an array of clustering results.

- k:

  number of clusters

## Value

cluster assignments for the consensus class

## References

Strehl, A., & Ghosh, J. (2002). Cluster ensembles—a knowledge reuse
framework for combining multiple partitions. Journal of machine learning
research, 3(Dec), 583-617.

## See also

Other consensus functions:
[`LCA()`](https://alinetalhouk.github.io/diceR/reference/LCA.md),
[`LCE()`](https://alinetalhouk.github.io/diceR/reference/LCE.md),
[`k_modes()`](https://alinetalhouk.github.io/diceR/reference/k_modes.md),
[`majority_voting()`](https://alinetalhouk.github.io/diceR/reference/majority_voting.md)

## Author

Derek Chiu

## Examples

``` r
data(hgsc)
dat <- hgsc[1:100, 1:50]
x <- consensus_cluster(dat, nk = 4, reps = 4, algorithms = c("hc", "diana"),
progress = FALSE)
CSPA(x, k = 4)
#>   [1] 1 1 2 2 2 1 1 2 2 1 2 2 1 1 3 2 1 1 2 1 1 1 1 1 2 1 1 1 1 1 1 1 1 2 1 1 2
#>  [38] 2 3 2 1 2 2 2 2 2 2 2 1 1 2 2 1 1 1 2 2 1 2 1 1 1 2 1 1 4 1 2 2 1 2 1 1 2
#>  [75] 1 2 1 1 3 1 2 1 2 2 1 2 1 1 3 2 2 2 1 2 2 1 1 2 1 1
```
