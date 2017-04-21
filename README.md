
<!-- README.md is generated from README.Rmd. Please edit that file -->
diceR
=====

[![Travis-CI Build Status](https://travis-ci.org/AlineTalhouk/diceR.svg?branch=master)](https://travis-ci.org/AlineTalhouk/diceR) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/AlineTalhouk/diceR?branch=master&svg=true)](https://ci.appveyor.com/project/dchiu911/diceR) [![Coverage Status](https://codecov.io/gh/AlineTalhouk/diceR/branch/master/graph/badge.svg)](https://codecov.io/gh/AlineTalhouk/diceR?branch=master)

The goal of `diceR` is to provide pipelines for generating diverse cluster ensembles in R.

Installation
------------

You can install `diceR` from github with:

``` r
# install.packages("devtools")
devtools::install_github("AlineTalhouk/diceR")
```

Example
-------

This basic example shows you how to use the main function of the package, `dice()`. A data matrix `hgsc` is partitioned into (a range of) `nk` clusters over `reps` bootstrap subsamples using each of the clustering `algorithms`. Clustering assignments are aggregated by the `cons.funs`.

``` r
library(diceR)
data(hgsc)
obj <- dice(hgsc, nk = 4, reps = 5, algorithms = c("hc", "diana"),
            cons.funs = c("kmodes", "majority"))
```

The first few cluster assignments are shown below:

``` r
knitr::kable(head(obj$clusters))
```

|                      |  kmodes|  majority|
|----------------------|-------:|---------:|
| TCGA.04.1331\_PRO.C5 |       3|         3|
| TCGA.04.1332\_MES.C1 |       3|         3|
| TCGA.04.1336\_DIF.C4 |       1|         3|
| TCGA.04.1337\_MES.C1 |       3|         3|
| TCGA.04.1338\_MES.C1 |       3|         3|
| TCGA.04.1341\_PRO.C5 |       3|         3|

You can also compare the base `algorithms` with the `cons.funs` using internal evaluation indices:

``` r
knitr::kable(obj$indices$internal$`4`)
```

| Algorithms       |   c\_index|  calinski\_harabasz|  davies\_bouldin|       dunn|  mcclain\_rao|       pbm|    sd\_dis|  ray\_turi|        tau|      gamma|    g\_plus|  silhouette|     s\_dbw|  Compactness|  Connectivity|
|:-----------------|----------:|-------------------:|----------------:|----------:|-------------:|---------:|----------:|----------:|----------:|----------:|----------:|-----------:|----------:|------------:|-------------:|
| HC\_Euclidean    |  0.3122823|            4.945499|         3.100302|  0.3025234|     0.8237540|  38.34704|  0.1795670|  3.0886000|  0.1992999|  0.5598731|  0.0278858|   0.0300838|        NaN|     24.81662|      49.69405|
| DIANA\_Euclidean |  0.1639431|           51.332198|         3.037874|  0.3348103|     0.8077658|  32.92726|  0.2034291|  3.1687896|  0.4271483|  0.6216897|  0.0892952|   0.0700862|        NaN|     22.05147|     227.34841|
| kmodes           |  0.2020221|           39.127460|         1.563373|  0.3352598|     0.8254116|  49.27019|  0.1046540|  1.1356906|  0.3907289|  0.5528538|  0.1116735|         NaN|  0.7207352|     22.66419|     148.61865|
| majority         |  0.2458043|            5.645220|         1.379460|  0.4315581|     0.7781939|  96.93674|  0.0948754|  0.8261741|  0.2221915|  0.7330421|  0.0122634|         NaN|  0.7224928|     24.70600|      24.35079|

Pipeline
--------

This figure is a visual schematic of the pipeline that `dice()` implements.

![Ensemble Clustering pipeline.](man/figures/pipeline.png)
