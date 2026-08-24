# Impute missing values

Impute missing values from bootstrapped subsampling

## Usage

``` r
impute_missing(E, data, nk)
```

## Arguments

- E:

  4D array of clusterings from `consensus_cluster`. The number of rows
  is equal to the number of cases to be clustered, number of columns is
  equal to the clusterings obtained by different resamplings of the
  data, the third dimension are the different algorithms and the fourth
  dimension are cluster sizes.

- data:

  data matrix with samples as rows and genes/features as columns

- nk:

  cluster size to extract data for (single value)

## Value

If flattened matrix consists of more than one repetition, i.e. it isn't
a column vector, then the function returns a matrix of clusterings with
complete cases imputed using majority voting, and relabelled, for chosen
`k`.

## Details

The default output from `consensus_cluster` will undoubtedly contain
`NA` entries because each replicate chooses a random subset (with
replacement) of all samples. Missing values should first be imputed
using
[`impute_knn()`](https://alinetalhouk.github.io/diceR/reference/impute_knn.md).
Not all missing values are guaranteed to be imputed by KNN. See
[`class::knn()`](https://rdrr.io/pkg/class/man/knn.html) for details.
Thus, any remaining missing values are imputed using majority voting.

## See also

Other imputation functions:
[`impute_knn()`](https://alinetalhouk.github.io/diceR/reference/impute_knn.md)

## Author

Aline Talhouk

## Examples

``` r
data(hgsc)
dat <- hgsc[1:100, 1:50]
E <- consensus_cluster(dat, nk = 3:4, reps = 10, algorithms = c("hc", "km",
"sc"), progress = FALSE)
sum(is.na(E))
#> [1] 1200
E_imputed <- impute_missing(E, dat, 4)
sum(is.na(E_imputed))
#> [1] 0
```
