# Simulate and select null distributions on empirical gene-gene correlations

Using a principal component constructed from the sample space, we
simulate null distributions with univariate Normal distributions using
`pcn_simulate`. Then a subset of these distributions is chosen using
`pcn_select`.

## Usage

``` r
pcn_simulate(data, n.sim = 50)

pcn_select(data.sim, cl, type = c("rep", "range"), int = 5)
```

## Arguments

- data:

  data matrix with rows as samples, columns as features

- n.sim:

  The number of simulated datasets to simulate

- data.sim:

  an object from `pcn_simulate`

- cl:

  vector of cluster memberships

- type:

  select either the representative dataset ("rep") or a range of
  datasets ("range")

- int:

  every `int` data sets from median-ranked `data.sim` are taken.
  Defaults to 5.

## Value

`pcn_simulate` returns a list of length `n.sim`. Each element is a
simulated matrix using this "Principal Component Normal" (pcn)
procedure.

`pcn_select` returns a list with elements

- `ranks`: When `type = "range"`, ranks of each extracted dataset shown

- `ind`: index of representative simulation

- `dat`: simulation data representation of all in pcNormal

## Author

Derek Chiu

## Examples

``` r
set.seed(9)
A <- matrix(rnorm(300), nrow = 20)
pc.dat <- pcn_simulate(A, n.sim = 50)
cl <- sample(1:4, 20, replace = TRUE)
pc.select <- pcn_select(pc.dat, cl, "rep")
```
