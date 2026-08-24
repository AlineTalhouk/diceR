# Minimize Frobenius norm for between two matrices

Finds a permutation of a matrix such that its Frobenius norm with
another matrix is minimized.

## Usage

``` r
min_fnorm(A, B = diag(nrow(A)))
```

## Arguments

- A:

  data matrix we want to permute

- B:

  matrix whose distance with the permuted A we want to minimize. By
  default, `B <- diag(nrow(A))`, so the permutation maximizes the trace
  of A.

## Value

Permuted matrix such that it is the permutation of A closest to B

## Details

Finds the permutation P of A such that `||PA - B||` is minimum in
Frobenius norm. Uses the linear-sum assignment problem (LSAP) solver in
the package `clue`. The default B is the identity matrix of same
dimension, so that the permutation of A maximizes its trace. This
procedure is useful for constructing a confusion matrix when we don't
know the true class labels of a predicted class and want to compare to a
reference class.

## Author

Ravi Varadhan:
https://stat.ethz.ch/pipermail/r-help/2010-April/236664.html

## Examples

``` r
set.seed(1)
A <- matrix(sample(1:25, size = 25, rep = FALSE), 5, 5)
min_fnorm(A)
#> $pmat
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]   25   11   16    9    8
#> [2,]    1   22   19   17    3
#> [3,]    2    5   23   20   24
#> [4,]    4   14   10   15   13
#> [5,]    7   18    6   12   21
#> 
#> $perm
#> Optimal assignment:
#> 1 => 1, 2 => 4, 3 => 5, 4 => 2, 5 => 3
#> 
#> $ord
#> [1] 1 4 5 2 3
#> 
```
