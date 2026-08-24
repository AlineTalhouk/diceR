# K-Nearest Neighbours imputation

The non-missing cases indicate the training set, and missing cases
indicate the test set.

## Usage

``` r
impute_knn(x, data, seed = 123456)
```

## Arguments

- x:

  clustering object

- data:

  data matrix

- seed:

  random seed for knn imputation reproducibility

## Value

An object with (potentially not all) missing values imputed with
K-Nearest Neighbours.

## Note

We consider 5 nearest neighbours and the minimum vote for definite
decision is 3.

## See also

Other imputation functions:
[`impute_missing()`](https://alinetalhouk.github.io/diceR/reference/impute_missing.md)

## Author

Aline Talhouk

## Examples

``` r
data(hgsc)
dat <- hgsc[1:100, 1:50]
x <- consensus_cluster(dat, nk = 4, reps = 4, algorithms = c("km", "hc",
"diana"), progress = FALSE)
x <- apply(x, 2:4, impute_knn, data = dat, seed = 1)
```
