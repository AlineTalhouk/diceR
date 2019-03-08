# diceR (development version)

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
