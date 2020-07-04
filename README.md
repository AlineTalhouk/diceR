
<!-- README.md is generated from README.Rmd. Please edit that file -->

# diceR

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/AlineTalhouk/diceR.svg?branch=master)](https://travis-ci.org/AlineTalhouk/diceR)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/AlineTalhouk/diceR?branch=master&svg=true)](https://ci.appveyor.com/project/AlineTalhouk/diceR)
[![Codecov test
coverage](https://codecov.io/gh/AlineTalhouk/diceR/branch/master/graph/badge.svg)](https://codecov.io/gh/AlineTalhouk/diceR?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/diceR)](https://CRAN.R-project.org/package=diceR)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/diceR?color=orange)](https://r-pkg.org/pkg/diceR)
<!-- badges: end -->

## Overview

The goal of `diceR` is to provide a systematic framework for generating
diverse cluster ensembles in R. There are a lot of nuances in cluster
analysis to consider. We provide a process and a suite of functions and
tools to implement a systematic framework for cluster discovery, guiding
the user through the generation of a diverse clustering solutions from
data, ensemble formation, algorithm selection and the arrival at a final
consensus solution. We have additionally developed visual and analytical
validation tools to help with the assessment of the final result. We
implemented a wrapper function `dice()` that allows the user to easily
obtain results and assess them. Thus, the package is accessible to both
end user with limited statistical knowledge. Full access to the package
is available for informaticians and statisticians and the functions are
easily expanded. More details can be found in our companion paper
published at [BMC
Bioinformatics](https://doi.org/10.1186/s12859-017-1996-y).

## Installation

You can install `diceR` from CRAN with:

``` r
install.packages("diceR")
```

Or get the latest development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("AlineTalhouk/diceR")
```

## Example

The following example shows how to use the main function of the package,
`dice()`. A data matrix `hgsc` contains a subset of gene expression
measurements of High Grade Serous Carcinoma Ovarian cancer patients from
the Cancer Genome Atlas publicly available datasets. Samples as rows,
features as columns. The function below runs the package through the
`dice()` function. We specify (a range of) `nk` clusters over `reps`
subsamples of the data containing 80% of the full samples. We also
specify the clustering `algorithms` to be used and the ensemble
functions used to aggregated them in `cons.funs`.

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

|                      | kmodes | majority |
| :------------------- | -----: | -------: |
| TCGA.04.1331\_PRO.C5 |      3 |        3 |
| TCGA.04.1332\_MES.C1 |      3 |        3 |
| TCGA.04.1336\_DIF.C4 |      1 |        3 |
| TCGA.04.1337\_MES.C1 |      3 |        3 |
| TCGA.04.1338\_MES.C1 |      3 |        3 |
| TCGA.04.1341\_PRO.C5 |      3 |        3 |

You can also compare the base `algorithms` with the `cons.funs` using
internal evaluation indices:

``` r
knitr::kable(obj$indices$ii$`4`)
```

| Algorithms       | calinski\_harabasz |      dunn |      pbm |       tau |     gamma |  c\_index | davies\_bouldin | mcclain\_rao |   sd\_dis | ray\_turi |   g\_plus | silhouette | s\_dbw | Compactness | Connectivity |
| :--------------- | -----------------: | --------: | -------: | --------: | --------: | --------: | --------------: | -----------: | --------: | --------: | --------: | ---------: | -----: | ----------: | -----------: |
| HC\_Euclidean    |           4.945499 | 0.3025234 | 38.34704 | 0.1992999 | 0.5598731 | 0.3122823 |        2.585780 |    0.8237540 | 0.1795670 | 3.0886000 | 0.0278858 |  0.0300838 |    NaN |    24.81662 |     49.69405 |
| DIANA\_Euclidean |          51.332198 | 0.3348103 | 32.92726 | 0.4271483 | 0.6216897 | 0.1639431 |        2.782694 |    0.8077658 | 0.2034291 | 3.1687896 | 0.0892952 |  0.0700862 |    NaN |    22.05147 |    227.34841 |
| kmodes           |          39.127460 | 0.3352598 | 49.27019 | 0.3907289 | 0.5528538 | 0.2020221 |        1.669753 |    0.8254116 | 0.1046540 | 1.1356906 | 0.1116735 |        NaN |    NaN |    22.66419 |    148.61865 |
| majority         |           5.645220 | 0.4315581 | 96.93674 | 0.2221915 | 0.7330421 | 0.2458043 |        1.393811 |    0.7781939 | 0.0948754 | 0.8261741 | 0.0122634 |        NaN |    NaN |    24.70600 |     24.35079 |

## Pipeline

This figure is a visual schematic of the pipeline that `dice()`
implements.

![Ensemble Clustering pipeline.](man/figures/pipeline.png)

Please visit the
[overview](https://alinetalhouk.github.io/diceR/articles/overview.html "diceR overview")
page for more detail.
