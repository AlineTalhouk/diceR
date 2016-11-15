
<!-- README.md is generated from README.Rmd. Please edit that file -->
diceR
=====

[![Travis-CI Build Status](https://travis-ci.com/AlineTalhouk/diceR.svg?token=R9saDTyWyg3zsFXfoa1H&branch=master)](https://travis-ci.com/AlineTalhouk/diceR) [![Coverage Status](https://codecov.io/gh/AlineTalhouk/diceR/branch/master/graph/badge.svg)](https://codecov.io/gh/AlineTalhouk/diceR?branch=master)

The goal of `diceR` is to provide pipelines for generating diverse cluster ensembles in R.

Installation
------------

You can install `diceR` from github with:

``` r
# install.packages("devtools")
# devtools::install_github("AlineTalhouk/diceR")
```

Example
-------

This basic example shows you how to use the main function of the package, `dice()`. A data matrix `dat` is partitioned into (a range of) `nk` clusters over `reps` bootstrap subsamples using each of the clustering `algorithms`. Clustering assignments are aggregated by the `cons.funs`.

``` r
library(diceR)
data(hgsc)
dat <- t(hgsc[, -1])
obj <- dice(dat, nk = 4, reps = 5, algorithms = c("hcAEucl", "hcDianaEucl"),
            cons.funs = c("kmodes", "majority"))
```

The first few cluster assignments are shown below:

``` r
knitr::kable(head(obj$clusters))
```

|  kmodes|  majority|
|-------:|---------:|
|       4|         1|
|       4|         1|
|       1|         1|
|       4|         1|
|       4|         1|
|       1|         1|

You can also compare the base `algorithms` with the `cons.funs` using internal evaluation indices:

``` r
knitr::kable(obj$indices$internal$`4`)
```

| Algorithms  |   c\_index|  calinski\_harabasz|  davies\_bouldin|       dunn|  mcclain\_rao|        pbm|    sd\_dis|  ray\_turi|        tau|      gamma|    g\_plus|  Compactness|  Connectivity|
|:------------|----------:|-------------------:|----------------:|----------:|-------------:|----------:|----------:|----------:|----------:|----------:|----------:|------------:|-------------:|
| hcAEucl     |  0.2108494|            19.41823|         3.292021|  0.3107600|     0.8077705|   20.76198|  0.1835442|   4.136054|  0.3829042|  0.6339045|  0.0667878|     23.84512|      90.57857|
| hcDianaEucl |  0.2197691|            52.63775|         2.201666|  0.3314038|     0.8402420|   40.85617|  0.1382588|   1.625841|  0.3374781|  0.4982350|  0.1151047|     21.96780|     245.91349|
| kmodes      |  0.1535169|            65.35652|         1.984004|  0.3403104|     0.8019908|   42.90585|  0.1771115|   2.148789|  0.4044636|  0.6382378|  0.0726419|     21.39036|     269.61310|
| majority    |  0.2272101|            27.50499|         2.038839|  0.3096773|     0.8148563|  129.67108|  0.1138123|   1.060246|  0.3620268|  0.6068949|  0.0699412|     23.84113|      72.42857|

Pipeline
--------

This figure is a visual schematic of the pipeline that `dice()` implements.

![Caption for the picture.](inst/img/pipeline.png)
