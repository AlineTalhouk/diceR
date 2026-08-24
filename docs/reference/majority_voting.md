# Majority voting

Combine clustering results using majority voting.

## Usage

``` r
majority_voting(E, is.relabelled = TRUE)
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

## Value

a vector of cluster assignments based on majority voting

## Details

Combine clustering results generated using different algorithms and
different data perturbations by majority voting. The class of a sample
is the cluster label which was selected most often across algorithms and
subsamples.

## See also

Other consensus functions:
[`CSPA()`](https://alinetalhouk.github.io/diceR/reference/CSPA.md),
[`LCA()`](https://alinetalhouk.github.io/diceR/reference/LCA.md),
[`LCE()`](https://alinetalhouk.github.io/diceR/reference/LCE.md),
[`k_modes()`](https://alinetalhouk.github.io/diceR/reference/k_modes.md)

## Author

Aline Talhouk

## Examples

``` r
data(hgsc)
dat <- hgsc[1:100, 1:50]
cc <- consensus_cluster(dat, nk = 4, reps = 6, algorithms = "pam", progress =
FALSE)
table(majority_voting(cc[, , 1, 1, drop = FALSE], is.relabelled = FALSE))
#> 
#>  1  2  3  4 
#> 40 32  9 19 
```
