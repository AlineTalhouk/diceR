# Consensus matrix

Returns the (weighted) consensus matrix given a data matrix

## Usage

``` r
consensus_matrix(data, weights = NULL)
```

## Arguments

- data:

  data matrix has rows as samples, columns as replicates

- weights:

  a vector of weights for each algorithm used in meta-consensus
  clustering. Must have `length(weights)` equal to `ncol(data)`.

## Value

a consensus matrix

## Details

Given a vector of cluster assignments, we first calculate the
connectivity matrix and indicator matrix. A connectivity matrix has a 1
if both samples are in the same cluster, and 0 otherwise. An indicator
matrix has a 1 if both samples were selected to be used in a subsample
of a consensus clustering algorithm, and 0 otherwise. Summation of
connectivity matrices and indicator matrices is performed over different
subsamples of the data. The consensus matrix is calculated by dividing
the aggregated connectivity matrices by the aggregated indicator
matrices.

If a meta-consensus matrix is desired, where consensus classes of
different clustering algorithms are aggregated, we can construct a
weighted meta-consensus matrix using `weights`.

## Note

When consensus is calculated over bootstrap samples, not every sample is
used in each replication. Thus, there will be scenarios where two
samples are never chosen together in any bootstrap samples. This
typically happens when the number of replications is small. The
coordinate in the consensus matrix for such pairs of samples is `NaN`
from a 0 / 0 computation. These entries are coerced to 0.

## Author

Derek Chiu

## Examples

``` r
set.seed(2)
x <- replicate(100, rbinom(100, 4, 0.2))
w <- rexp(100)
w <- w / sum(w)
cm1 <- consensus_matrix(x)
cm2 <- consensus_matrix(x, weights = w)
```
