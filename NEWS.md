# diceR 0.3.0.9000

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
