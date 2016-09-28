
<!-- README.md is generated from README.Rmd. Please edit that file -->
ClusterDuck
===========

The goal of `ClusterDuck` is to provide a automated pipeline for ensemble clustering.

Installation
------------

You can install ClusterDuck from github with:

``` r
# install.packages("devtools")
# devtools::install_github("AlineTalhouk/ClusterDuck")
```

Example
-------

This is a basic example which shows you how to use the main workhouse function of the package, `ConClust()`:

``` r
library(ClusterDuck)
data(hgsc)
x <- ConClust(hgsc, k = 4, reps = 10, method = "hcAEucl", save = FALSE)
#> 
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |======                                                           |  10%
  |                                                                       
  |=============                                                    |  20%
  |                                                                       
  |====================                                             |  30%
  |                                                                       
  |==========================                                       |  40%
  |                                                                       
  |================================                                 |  50%
  |                                                                       
  |=======================================                          |  60%
  |                                                                       
  |==============================================                   |  70%
  |                                                                       
  |====================================================             |  80%
  |                                                                       
  |==========================================================       |  90%
  |                                                                       
  |=================================================================| 100%
head(x[, , 1])
#>                     R1 R2 R3 R4 R5 R6 R7 R8 R9 R10
#> TCGA.04.1331_PRO.C5 NA  1  1  1  1  1 NA  1  1   1
#> TCGA.04.1332_MES.C1 NA  1  1  1  1  1 NA  1  1   1
#> TCGA.04.1336_DIF.C4 NA  1  1  1  1  1 NA  1  1   1
#> TCGA.04.1337_MES.C1 NA  1  1  1  1  1  1 NA  1   1
#> TCGA.04.1338_MES.C1  1  1  1  1  1  1  1  1  1   1
#> TCGA.04.1341_PRO.C5  1  1  1  1 NA  1 NA  1  1   1
```
