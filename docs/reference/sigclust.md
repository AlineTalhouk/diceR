# Significant Testing of Clustering Results

Uses the SigClust K-Means algorithm to assess significance of clustering
results.

## Usage

``` r
sigclust(x, k, nsim, nrep = 1, labflag = 0, label = 0, icovest = 2)
```

## Arguments

- x:

  data matrix, samples are rows and features are columns

- k:

  cluster size to test against

- nsim:

  number of simulations

- nrep:

  See
  [`sigclust::sigclust()`](https://rdrr.io/pkg/sigclust/man/sigclust.html)
  for details.

- labflag:

  See
  [`sigclust::sigclust()`](https://rdrr.io/pkg/sigclust/man/sigclust.html)
  for details.

- label:

  true class label. See
  [`sigclust::sigclust()`](https://rdrr.io/pkg/sigclust/man/sigclust.html)
  for details.

- icovest:

  type of covariance matrix estimation

## Value

An object of class `sigclust`. See
[`sigclust::sigclust()`](https://rdrr.io/pkg/sigclust/man/sigclust.html)
for details.

## Details

This function is a wrapper for the original
[`sigclust::sigclust()`](https://rdrr.io/pkg/sigclust/man/sigclust.html),
except that an additional parameter `k` is allows testing against any
number of clusters. In addition, the default type of covariance
estimation is also different.

## References

Liu, Yufeng, Hayes, David Neil, Nobel, Andrew and Marron, J. S, 2008,
*Statistical Significance of Clustering for High-Dimension, Low-Sample
Size Data*, *Journal of the American Statistical Association*
**103**(483) 1281–1293.

## Author

Hanwen Huang: <hanwenh@email.unc.edu>; Yufeng Liu:
<yfliu@email.unc.edu>; J. S. Marron: <marron@email.unc.edu>

## Examples

``` r
data(hgsc)
dat <- hgsc[1:100, 1:50]
nk <- 4
cc <- consensus_cluster(dat, nk = nk, reps = 5, algorithms = "pam",
progress = FALSE)
cl.mat <- consensus_combine(cc, element = "class")
lab <- cl.mat$`4`[, 1]
set.seed(1)
str(sigclust(x = dat, k = nk, nsim = 50, labflag = 1, label = lab))
#> Formal class 'sigclust' [package "sigclust"] with 10 slots
#>   ..@ raw.data  : num [1:100, 1:50] -0.0107 -0.7107 0.8815 -1.0851 -0.9322 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:100] "TCGA.04.1331_PRO.C5" "TCGA.04.1332_MES.C1" "TCGA.04.1336_DIF.C4" "TCGA.04.1337_MES.C1" ...
#>   .. .. ..$ : chr [1:50] "ABAT" "ABHD2" "ACTB" "ACTR2" ...
#>   ..@ veigval   : num [1:50] 11.81 4.51 2.66 2.29 1.84 ...
#>   ..@ vsimeigval: num [1:50] 11.81 4.51 2.66 2.29 1.84 ...
#>   ..@ simbackvar: num 0.42
#>   ..@ icovest   : num 2
#>   ..@ nsim      : num 50
#>   ..@ simcindex : num [1:50] 0.647 0.673 0.65 0.652 0.59 ...
#>   ..@ pval      : num 0.76
#>   ..@ pvalnorm  : num 0.776
#>   ..@ xcindex   : num 0.667
```
