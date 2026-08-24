# Latent Class Analysis

Combine clustering results using latent class analysis.

## Usage

``` r
LCA(E, is.relabelled = TRUE, seed = 1)
```

## Arguments

- E:

  a matrix of clusterings with number of rows equal to the number of
  cases to be clustered, number of columns equal to the clustering
  obtained by different resampling of the data, and the third dimension
  are the different algorithms. Matrix may already be two-dimensional.

- is.relabelled:

  logical; if `FALSE` the data will be relabelled using the first
  clustering as the reference.

- seed:

  random seed for reproducibility

## Value

a vector of cluster assignments based on LCA

## See also

Other consensus functions:
[`CSPA()`](https://alinetalhouk.github.io/diceR/reference/CSPA.md),
[`LCE()`](https://alinetalhouk.github.io/diceR/reference/LCE.md),
[`k_modes()`](https://alinetalhouk.github.io/diceR/reference/k_modes.md),
[`majority_voting()`](https://alinetalhouk.github.io/diceR/reference/majority_voting.md)

## Author

Derek Chiu

## Examples

``` r
data(hgsc)
dat <- hgsc[1:100, 1:50]
cc <- consensus_cluster(dat, nk = 4, reps = 6, algorithms = "pam", progress =
FALSE)
table(LCA(cc[, , 1, 1, drop = FALSE], is.relabelled = FALSE))
#> 
#>  1  2  3  4 
#> 21 38  9 32 
```
