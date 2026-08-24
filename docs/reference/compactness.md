# Compactness Measure

Compute the compactness validity index for a clustering result.

## Usage

``` r
compactness(data, labels)
```

## Arguments

- data:

  a dataset with rows as observations, columns as variables

- labels:

  a vector of cluster labels from a clustering result

## Value

the compactness score

## Details

This index is agnostic to any reference clustering results, calculating
cluster performance on the basis of compactness and separability.
Smaller values indicate a better clustering structure.

## References

MATLAB function `valid_compactness` by Simon Garrett in LinkCluE

## Author

Derek Chiu

## Examples

``` r
set.seed(1)
E <- matrix(rep(sample(1:4, 1000, replace = TRUE)), nrow = 100, byrow =
              FALSE)
set.seed(1)
dat <- as.data.frame(matrix(runif(1000, -10, 10), nrow = 100, byrow = FALSE))
compactness(dat, E[, 1])
#> [1] 25.32475
```
