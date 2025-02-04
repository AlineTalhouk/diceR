# diceR 3.0.0

* There were errors running tests on clang-UBSAN and gcc-UBSAN builds, due to usage of functions from `clusterSim`. Reinstate internal and external functions from `clusterCrit` package, which is no longer archived. 
* Use new parameter `verbose` in `dice()` to control console printouts of main tasks being performed instead of using `progress`
* Refactor `graph_heatmap()` to use `pheatmap::pheatmap()` since `NMF::aheatmap()` throws a `gridPLT()` error whenever it is run in a script or interactively, but not in R markdown documents. Thus there are run time errors when used in unit tests and examples. See renozao/NMF#65
* Remove deprecated function usage

# diceR 2.2.0

* New `abs` argument in `consensus_cluster()`: control whether to apply absolute value to Spearman and Pearson correlation matrices before subtracting from one (@tiagochst, #161)
* New distance matrix option, `distance = "pearson"`
* Use `blockcluster` instead of `mixedClust` as the latter is now archived
* Update tests in `consensus_cluster()` that skip when suggested packages are not installed

# diceR 2.1.0

* Sort cluster sizes `k` correctly in relative change in area under CDF curve (@IgnatiusPang, #167) and consensus matrix CDF graphs
* Replace deprecated `aes_` calls with tidyeval idioms
* Pass `p.item` and `seed.data` arguments to `dice()` (#162, #165)

# diceR 2.0.0

Internal and external validity indices were refactored to avoid using helper functions from the `clusterCrit` package, which is scheduled to be archived. Please adapt your code if it extracts deprecated validity indices, as described below.

* Calinski-Harabasz index now calculated using `clusterSim::index.G1()`
* Dunn index now calculated using `clValid::dunn()`
* Gamma index now calculated using `clusterSim::index.G2()`
* C-index now calculated using `clusterSim::index.C()`
* Davies-Bouldin now calculated using `clusterSim::index.DB()`
* SD index now calculated _correctly_ using `clv::clv.SD()` and helper functions from `clv`. Previously only the total separation between clusters was returned.
* S_Dbw index now calculated using `clv::clv.SDbw()` and helper functions from `clv`
* Rousseeuw's Silhouette now calculated using `clusterSim::index.S()`
* PBM, Tau, McClain-Rao, Ray-Turi, and G-plus indices were removed as equivalent implementations from other packages were not found. They may be reinstated in the future.
* All external validity indices are now calculated manually using counts from the concordance matrix (Hubert, Jaccard, McNemar, and Rand indices)
* Best index value (maximum or minimum) now calculated manually 

# diceR 1.2.2

* Pass `lower` and `upper` parameters from `PAC()` to `consensus_evaluate()` (#160)

# diceR 1.2.1

* References added for `k_modes()` and `CSPA()` consensus functions (#157)
* Update roxygen2 to avoid HTML5 documentation notes
* Use `mixedClust` instead of `blockcluster` for co-clustering since the latter keeps getting archived

# diceR 1.2.0

* Use `testthat::skip_if_not_installed()` to run tests conditionally when using packages in Suggests
* Use roxygen tag `@examplesIf` and `rlang::is_installed()` to run examples conditionally when using packages in Suggests

# diceR 1.1.0

* Reinstate `blockcluster` package as it is under active maintenance again

# diceR 1.0.4

* Suppress new names messages from transformed NMF data

* Flattened matrices include 4th dimension of clustering array

# diceR 1.0.3

* Add package logo using `hexSticker`

* Package `blockcluster` has been archived, remove from Suggests

# diceR 1.0.2

* Remove deprecated `context()` in tests

* Add `RColorBrewer` to Suggests because it is cross referenced in `?graphs`

* Add trailing slash for URLs in DESCRIPTION

* Remove `tibble` from Imports, no longer used

# diceR 1.0.1

* Suppress warnings when `clValid::connectivity()` is called regarding comparison with more than one class. Since R-4.0.0, a `matrix` object also inherits from class `array`

* In `algii_heatmap()`, the object `ii` already has row names passed from columns which are named vectors (issue also identified in #148, thanks @romainfrancois)

# diceR 1.0.0

## Decreased dependencies

The following steps were taken to minimize dependencies and ensure that `diceR` can still run on R 3.5:

* Removed `cli` and `RColorBrewer` from Imports

* Moved `apcluster`, `blockcluster`, `cluster`, `dbscan`, `e1071`, `kernlab`, and `kohonen` to `Suggests`, use their specific clustering algorithms conditionally. `mclust` needs to be in `Imports` because `mclust::mclustBIC()` needs to be imported

* Moved `sigclust` to `Suggests`, use within `sigclust()` conditionally

* Moved `progress` to `Suggests`, use within `consensus_cluster()` conditionally

* Moved `poLCA` to `Suggests`, use within `dice()` conditionally

* Moved `Rtsne` to `Suggests`, use within `prepare_data()` conditionally

* Removed old dependency `grDevices` from `Imports`

* Set minimum version to R (>= 3.5) for `klaR` dependency `questionr`

* In `ev_confmat()`, use `yardstick::conf_mat()` instead of `caret::confusionMatrix()`. `caret` has many dependencies, so best to avoid using it

* In `graph_heatmap()`, use `NMF::aheatmap()` instead of `gplots::heatmap.2()`. `gplots` depends on `caTools`, which now relies on R (>= 3.6)

* In `consensus_cluster()`, use `stringr::str_to_title()` instead of `Hmisc::capitalize()`. `Hmisc` depends on `latticeExtra`, which now relies on R (>= 3.6)

* In `graph_delta_area()`, use base solution instead of `flux::auc()`. `flux` also depends on `caTools`

* In `prepare_data()`, use own implementation of `quantable::robustscale()` with all of the former function's defaults. `quantable` also depends on `caTools`

* Specify Bioconductor installation on Travis and AppVeyor since `NMF` now Imports `Biobase`

## Minor improvements and bug fixes

* Remove `suppressWarnings(RNGversion("3.5.0"))` after updating R version

* Run `LCA()` unit test on imputed clustering object

* Remove internal validity measures with any `Inf` entries for `consensus_reweigh()`

* Use a cleaner, more robust method of removing `Rplots.pdf` after running `test-graphs.R`

* Ensure column binding with `purrr::map_dfc()` in `consensus_rank()`

* Replaced `dplyr::bind_cols()` with `purrr::flatten_dfc()` to suppress warning "Outer names are only allowed for unnamed scalar" in `get_cdf()`

* update roxygen and docs

# diceR 0.6.0

* Remove deprecated `dplyr` functions and use `.data` pronoun

* k-means clustering should not support distance matrices as input (@jerryji1993, #139)

* Add LCA as a consensus function (@philstraforelli, #137)

# diceR 0.5.2

* Fix `length > 1 in coercion to logical` error in `consensus_evaluate()` due to comparisons using `||` operator

* Add `suppressWarnings(RNGversion("3.5.0"))` before call to `set.seed()` in examples, tests, and vignette to use old RNG sampling

* Use `.covrignore` to exclude `zzz.R` from being considered in code coverage

* Use `dplyr` version >= 0.7.5 to ensure `bind_rows()` works

* Fixed bug where scaled matrix using the "robust" method in `prepare_data()` was nested in the `data` element (@AlineTalhouk, #134)

# diceR 0.5.1

* Add parameter `hc.method` in `dice` and `consensus_cluster` to pass to `method` parameter in `stats::hclust` (@JakeNel28, #130)

* Remove dependencies on `largeVis`: package will be archived

# diceR 0.5.0

* Revert back to using `NMF` since `NNLM` has been archived and `NMF` is back in active maintenance.

* Choose fuzzifier m in `cmeans` using Equation 5 from https://academic.oup.com/bioinformatics/article/26/22/2841/227572 (thanks @Asduveneck)

# diceR 0.4.0

* Replace all code that depended on `NMF` with `NNLM` and `pheatmap`: CRAN notified that `NMF` will be archived because of inactive maintenance

* Update `.yml` files default templates

# diceR 0.3.2

* Fix bug in `consensus_cluster()` when custom algorithms were excluded from output (thanks @phiala)

* Use markdown language for documentation

* Various performance improvements and code simplifications

# diceR 0.3.1

* Suppress success/fail message printout and fix input data to be matrix for block clustering

* Fix bug in `algii_heatmap()` when `k.method = "all"` in `dice()`

* Fix bug in calculating internal indices when data has categorical variables (thanks Kurt Salmela)

# diceR 0.3.0

* Updated object output names in `consensus_evaluate()`

* Fix unit test in `test-dice.R` for R-devel

* Add internal function: ranked algorithms vs internal validity indices heatmap graph

* Fix bugs in `graph_cdf()`, `graph_tracking()` when only one k selected

* Progress messages in `dice()`

* Fix bug in `consensus_evaluate()` when algorithm has `NA` for all PAC values

# diceR 0.2.0

* New dimension reduction methods: t-SNE, largeVis (@dustin21)

* Better annotated progress bar using `progress` package

* Speed up the operation that transforms a matrix to become "NMF-ready"

* Simplify saving mechanism in `consensus_cluster()` such that only `file.name` needs to be specified, and the `save` parameter has been removed

* New algorithms: SOM, Fuzzy C-Means, DBSCAN (@dustin21, #118)

* Added significance testing section to vignette

* Fixed direction of optimization: compactness and connectivity should be minimized

# diceR 0.1.0

* New submission to CRAN accepted on June 21, 2017
