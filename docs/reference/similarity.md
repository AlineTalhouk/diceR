# Similarity Matrices

`cts` computes the connected triple based similarity matrix, `srs`
computes the simrank based similarity matrix, and `asrs` computes the
approximated simrank based similarity matrix.

## Usage

``` r
cts(E, dc)

srs(E, dc, R)

asrs(E, dc)
```

## Arguments

- E:

  an N by M matrix of cluster ensembles

- dc:

  decay factor, ranges from 0 to 1 inclusive

- R:

  number of iterations for `srs`

## Value

an N by N CTS, SRS, or ASRS matrix

## References

MATLAB functions cts, srs, asrs in package LinkCluE by Simon Garrett

## Author

Johnson Liu, Derek Chiu

## Examples

``` r
set.seed(1)
E <- matrix(rep(sample(1:4, 800, replace = TRUE)), nrow = 100)
CTS <- cts(E = E, dc = 0.8)
SRS <- srs(E = E, dc = 0.8, R = 3)
ASRS <- asrs(E = E, dc = 0.8)
purrr::walk(list(CTS, SRS, ASRS), str)
#>  num [1:100, 1:100] 1 0.888 0.792 0.811 0.715 ...
#>  num [1:100, 1:100] 1 0.06901 0.03048 0.04348 0.00553 ...
#>  num [1:100, 1:100] 1 0.638 0.634 0.638 0.653 ...
```
