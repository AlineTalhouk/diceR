# Prepare data for consensus clustering

Perform feature selection or dimension reduction to remove noise
variables.

## Usage

``` r
prepare_data(
  data,
  scale = TRUE,
  type = c("conventional", "robust", "tsne"),
  min.var = 1
)
```

## Arguments

- data:

  data matrix with rows as samples and columns as variables

- scale:

  logical; should the data be centered and scaled?

- type:

  if we use "conventional" measures (default), then the mean and
  standard deviation are used for centering and scaling, respectively.
  If "robust" measures are specified, the median and median absolute
  deviation (MAD) are used. Alternatively, we can apply "tsne" for
  dimension reduction.

- min.var:

  minimum variability measure threshold used to filter the feature space
  for only highly variable features. Only features with a minimum
  variability measure across all samples greater than `min.var` will be
  used. If `type = "conventional"`, the standard deviation is the
  measure used, and if `type = "robust"`, the MAD is the measure used.

## Value

dataset prepared for usage in `consensus_cluster`

## Details

We can apply a basic filtering method of feature selection that removes
variables with low signal and (optionally) scales before consensus
clustering. Or, we can use t-SNE dimension reduction to transform the
data to just two variables. This lower-dimensional embedding allows
algorithms such as hierarchical clustering to achieve greater
performance.

## Author

Derek Chiu

## Examples

``` r
set.seed(2)
x <- replicate(10, rnorm(100))
x.prep <- prepare_data(x)
dim(x)
#> [1] 100  10
dim(x.prep)
#> [1] 100   4
```
