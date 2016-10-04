
<!-- README.md is generated from README.Rmd. Please edit that file -->
diceR
=====

[![Travis-CI Build Status](https://travis-ci.com/AlineTalhouk/diceR.svg?branch=master)](https://travis-ci.com/AlineTalhouk/diceR) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/AlineTalhouk/diceR?branch=master&svg=true)](https://ci.appveyor.com/project/dchiu911/diceR) [![Coverage Status](https://img.shields.io/codecov/c/github/AlineTalhouk/diceR/master.svg)](https://codecov.io/github/AlineTalhouk/diceR?branch=master)

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

This is a basic example which shows you how to use the main workhouse function of the package, `ConClust()`:

``` r
library(diceR)
data(hgsc)
x <- ConClust(hgsc, k = 4, reps = 10, method = "hcAEucl", save = FALSE)
```

``` r
head(x[, , 1])
#>                     R1 R2 R3 R4 R5 R6 R7 R8 R9 R10
#> TCGA.04.1331_PRO.C5 NA  1  1  1  1  1 NA  1  1   1
#> TCGA.04.1332_MES.C1 NA  1  1  1  1  1 NA  1  1   1
#> TCGA.04.1336_DIF.C4 NA  1  1  1  1  1 NA  1  1   1
#> TCGA.04.1337_MES.C1 NA  1  1  1  1  1  1 NA  1   1
#> TCGA.04.1338_MES.C1  1  1  1  1  1  1  1  1  1   1
#> TCGA.04.1341_PRO.C5  1  1  1  1 NA  1 NA  1  1   1
```
