# Relabel classes to a standard

Relabel clustering categories to match to a standard by minimizing the
Frobenius norm between the two labels.

## Usage

``` r
relabel_class(pred.cl, ref.cl)
```

## Arguments

- pred.cl:

  vector of predicted cluster assignments

- ref.cl:

  vector of reference labels to match to

## Value

A vector of relabeled cluster assignments

## Author

Aline Talhouk

## Examples

``` r
set.seed(2)
pred <- sample(1:4, 100, replace = TRUE)
true <- sample(1:4, 100, replace = TRUE)
relabel_class(pred, true)
#>   [1] 3 2 1 1 4 4 3 3 3 4 3 3 1 2 3 2 1 4 4 1 2 2 2 4 2 3 1 3 4 4 2 1 3 1 4 4 3
#>  [38] 1 4 1 2 3 2 4 3 3 1 1 1 2 1 2 4 2 3 2 4 3 3 2 3 1 4 4 3 4 1 3 4 1 3 1 3 2
#>  [75] 1 1 2 4 4 4 4 3 1 2 1 1 1 4 2 2 3 2 4 1 3 4 3 3 3 1
```
