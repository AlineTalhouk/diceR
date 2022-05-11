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
