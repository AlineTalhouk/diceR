# Linkage Clustering Ensemble

Generate a cluster assignment from a CTS, SRS, or ASRS similarity
matrix.

## Usage

``` r
LCE(E, k, dc = 0.8, R = 10, sim.mat = c("cts", "srs", "asrs"))
```

## Arguments

- E:

  is an array of clustering results. An error is thrown if there are
  missing values.
  [`impute_missing()`](https://alinetalhouk.github.io/diceR/reference/impute_missing.md)
  can be used beforehand.

- k:

  requested number of clusters

- dc:

  decay constant for CTS, SRS, or ASRS matrix

- R:

  number of repetitions for SRS matrix

- sim.mat:

  similarity matrix; choices are "cts", "srs", "asrs".

## Value

a vector containing the cluster assignment from either the CTS, SRS, or
ASRS similarity matrices

## See also

Other consensus functions:
[`CSPA()`](https://alinetalhouk.github.io/diceR/reference/CSPA.md),
[`LCA()`](https://alinetalhouk.github.io/diceR/reference/LCA.md),
[`k_modes()`](https://alinetalhouk.github.io/diceR/reference/k_modes.md),
[`majority_voting()`](https://alinetalhouk.github.io/diceR/reference/majority_voting.md)

## Author

Johnson Liu

## Examples

``` r
data(hgsc)
dat <- hgsc[1:100, 1:50]
x <- consensus_cluster(dat, nk = 4, reps = 4, algorithms = c("km", "hc"),
progress = FALSE)
if (FALSE) { # \dontrun{
LCE(E = x, k = 4, sim.mat = "asrs")
} # }

x <- apply(x, 2:4, impute_knn, data = dat, seed = 1)
x_imputed <- impute_missing(x, dat, nk = 4)
LCE(E = x_imputed, k = 4, sim.mat = "cts")
#>   [1] 1 1 1 1 1 1 1 1 1 2 1 1 2 1 3 1 1 2 1 1 2 1 2 1 1 1 2 1 2 1 1 1 1 1 2 2 1
#>  [38] 1 4 1 2 1 1 1 1 1 1 1 1 2 1 1 2 2 2 1 1 1 1 1 2 2 1 1 2 2 2 1 1 2 1 2 1 1
#>  [75] 1 1 1 1 4 2 1 1 1 1 1 1 1 1 3 1 1 1 2 1 1 2 2 1 2 1
```
