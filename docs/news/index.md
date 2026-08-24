# Changelog

## diceR 3.2.0

- Reinstate `mixedClust` for block clustering now that it is back on
  CRAN and `blockcluster` has been archived

## diceR 3.1.0

CRAN release: 2025-06-19

- Remove deprecated helper functions package `clv` currently pending
  archival
- Improve documentation for parameter `k.method` in
  `consensus_evaluation()`

## diceR 3.0.0

CRAN release: 2025-02-05

- There were errors running tests on clang-UBSAN and gcc-UBSAN builds,
  due to usage of functions from `clusterSim`. Reinstate internal and
  external functions from `clusterCrit` package, which is no longer
  archived.
- Use new parameter `verbose` in
  [`dice()`](https://alinetalhouk.github.io/diceR/reference/dice.md) to
  control console printouts of main tasks being performed instead of
  using `progress`
- Refactor
  [`graph_heatmap()`](https://alinetalhouk.github.io/diceR/reference/graphs.md)
  to use
  [`pheatmap::pheatmap()`](https://rdrr.io/pkg/pheatmap/man/pheatmap.html)
  since [`NMF::aheatmap()`](https://rdrr.io/pkg/NMF/man/aheatmap.html)
  throws a `gridPLT()` error whenever it is run in a script or
  interactively, but not in R markdown documents. Thus there are run
  time errors when used in unit tests and examples. See renozao/NMF#65
- Remove deprecated function usage

## diceR 2.2.0

CRAN release: 2024-01-22

- New `abs` argument in
  [`consensus_cluster()`](https://alinetalhouk.github.io/diceR/reference/consensus_cluster.md):
  control whether to apply absolute value to Spearman and Pearson
  correlation matrices before subtracting from one
  ([@tiagochst](https://github.com/tiagochst),
  [\#161](https://github.com/AlineTalhouk/diceR/issues/161))
- New distance matrix option, `distance = "pearson"`
- Use `blockcluster` instead of `mixedClust` as the latter is now
  archived
- Update tests in
  [`consensus_cluster()`](https://alinetalhouk.github.io/diceR/reference/consensus_cluster.md)
  that skip when suggested packages are not installed

## diceR 2.1.0

CRAN release: 2023-09-28

- Sort cluster sizes `k` correctly in relative change in area under CDF
  curve ([@IgnatiusPang](https://github.com/IgnatiusPang),
  [\#167](https://github.com/AlineTalhouk/diceR/issues/167)) and
  consensus matrix CDF graphs
- Replace deprecated `aes_` calls with tidyeval idioms
- Pass `p.item` and `seed.data` arguments to
  [`dice()`](https://alinetalhouk.github.io/diceR/reference/dice.md)
  ([\#162](https://github.com/AlineTalhouk/diceR/issues/162),
  [\#165](https://github.com/AlineTalhouk/diceR/issues/165))

## diceR 2.0.0

CRAN release: 2023-03-11

Internal and external validity indices were refactored to avoid using
helper functions from the `clusterCrit` package, which is scheduled to
be archived. Please adapt your code if it extracts deprecated validity
indices, as described below.

- Calinski-Harabasz index now calculated using
  [`clusterSim::index.G1()`](https://rdrr.io/pkg/clusterSim/man/index.G1.html)
- Dunn index now calculated using
  [`clValid::dunn()`](https://rdrr.io/pkg/clValid/man/dunn.html)
- Gamma index now calculated using
  [`clusterSim::index.G2()`](https://rdrr.io/pkg/clusterSim/man/index.G2.html)
- C-index now calculated using
  [`clusterSim::index.C()`](https://rdrr.io/pkg/clusterSim/man/index.C.html)
- Davies-Bouldin now calculated using
  [`clusterSim::index.DB()`](https://rdrr.io/pkg/clusterSim/man/index.DB.html)
- SD index now calculated *correctly* using `clv::clv.SD()` and helper
  functions from `clv`. Previously only the total separation between
  clusters was returned.
- S_Dbw index now calculated using `clv::clv.SDbw()` and helper
  functions from `clv`
- Rousseeuw’s Silhouette now calculated using
  [`clusterSim::index.S()`](https://rdrr.io/pkg/clusterSim/man/index.S.html)
- PBM, Tau, McClain-Rao, Ray-Turi, and G-plus indices were removed as
  equivalent implementations from other packages were not found. They
  may be reinstated in the future.
- All external validity indices are now calculated manually using counts
  from the concordance matrix (Hubert, Jaccard, McNemar, and Rand
  indices)
- Best index value (maximum or minimum) now calculated manually

## diceR 1.2.2

CRAN release: 2022-09-29

- Pass `lower` and `upper` parameters from
  [`PAC()`](https://alinetalhouk.github.io/diceR/reference/PAC.md) to
  [`consensus_evaluate()`](https://alinetalhouk.github.io/diceR/reference/consensus_evaluate.md)
  ([\#160](https://github.com/AlineTalhouk/diceR/issues/160))

## diceR 1.2.1

CRAN release: 2022-08-16

- References added for
  [`k_modes()`](https://alinetalhouk.github.io/diceR/reference/k_modes.md)
  and [`CSPA()`](https://alinetalhouk.github.io/diceR/reference/CSPA.md)
  consensus functions
  ([\#157](https://github.com/AlineTalhouk/diceR/issues/157))
- Update roxygen2 to avoid HTML5 documentation notes
- Use `mixedClust` instead of `blockcluster` for co-clustering since the
  latter keeps getting archived

## diceR 1.2.0

CRAN release: 2022-05-13

- Use
  [`testthat::skip_if_not_installed()`](https://testthat.r-lib.org/reference/skip.html)
  to run tests conditionally when using packages in Suggests
- Use roxygen tag `@examplesIf` and
  [`rlang::is_installed()`](https://rlang.r-lib.org/reference/is_installed.html)
  to run examples conditionally when using packages in Suggests

## diceR 1.1.0

CRAN release: 2021-07-23

- Reinstate `blockcluster` package as it is under active maintenance
  again

## diceR 1.0.4

CRAN release: 2021-06-04

- Suppress new names messages from transformed NMF data

- Flattened matrices include 4th dimension of clustering array

## diceR 1.0.3

CRAN release: 2021-04-17

- Add package logo using `hexSticker`

- Package `blockcluster` has been archived, remove from Suggests

## diceR 1.0.2

CRAN release: 2021-03-18

- Remove deprecated `context()` in tests

- Add `RColorBrewer` to Suggests because it is cross referenced in
  [`?graphs`](https://alinetalhouk.github.io/diceR/reference/graphs.md)

- Add trailing slash for URLs in DESCRIPTION

- Remove `tibble` from Imports, no longer used

## diceR 1.0.1

CRAN release: 2021-01-30

- Suppress warnings when
  [`clValid::connectivity()`](https://rdrr.io/pkg/clValid/man/connectivity.html)
  is called regarding comparison with more than one class. Since
  R-4.0.0, a `matrix` object also inherits from class `array`

- In `algii_heatmap()`, the object `ii` already has row names passed
  from columns which are named vectors (issue also identified in
  [\#148](https://github.com/AlineTalhouk/diceR/issues/148), thanks
  [@romainfrancois](https://github.com/romainfrancois))

## diceR 1.0.0

CRAN release: 2020-07-07

### Decreased dependencies

The following steps were taken to minimize dependencies and ensure that
`diceR` can still run on R 3.5:

- Removed `cli` and `RColorBrewer` from Imports

- Moved `apcluster`, `blockcluster`, `cluster`, `dbscan`, `e1071`,
  `kernlab`, and `kohonen` to `Suggests`, use their specific clustering
  algorithms conditionally. `mclust` needs to be in `Imports` because
  [`mclust::mclustBIC()`](https://mclust-org.github.io/mclust/reference/mclustBIC.html)
  needs to be imported

- Moved `sigclust` to `Suggests`, use within
  [`sigclust()`](https://alinetalhouk.github.io/diceR/reference/sigclust.md)
  conditionally

- Moved `progress` to `Suggests`, use within
  [`consensus_cluster()`](https://alinetalhouk.github.io/diceR/reference/consensus_cluster.md)
  conditionally

- Moved `poLCA` to `Suggests`, use within
  [`dice()`](https://alinetalhouk.github.io/diceR/reference/dice.md)
  conditionally

- Moved `Rtsne` to `Suggests`, use within
  [`prepare_data()`](https://alinetalhouk.github.io/diceR/reference/prepare_data.md)
  conditionally

- Removed old dependency `grDevices` from `Imports`

- Set minimum version to R (\>= 3.5) for `klaR` dependency `questionr`

- In
  [`ev_confmat()`](https://alinetalhouk.github.io/diceR/reference/external_validity.md),
  use
  [`yardstick::conf_mat()`](https://yardstick.tidymodels.org/reference/conf_mat.html)
  instead of
  [`caret::confusionMatrix()`](https://rdrr.io/pkg/caret/man/confusionMatrix.html).
  `caret` has many dependencies, so best to avoid using it

- In
  [`graph_heatmap()`](https://alinetalhouk.github.io/diceR/reference/graphs.md),
  use [`NMF::aheatmap()`](https://rdrr.io/pkg/NMF/man/aheatmap.html)
  instead of
  [`gplots::heatmap.2()`](https://talgalili.github.io/gplots/reference/heatmap.2.html).
  `gplots` depends on `caTools`, which now relies on R (\>= 3.6)

- In
  [`consensus_cluster()`](https://alinetalhouk.github.io/diceR/reference/consensus_cluster.md),
  use
  [`stringr::str_to_title()`](https://stringr.tidyverse.org/reference/case.html)
  instead of
  [`Hmisc::capitalize()`](https://rdrr.io/pkg/Hmisc/man/capitalize.html).
  `Hmisc` depends on `latticeExtra`, which now relies on R (\>= 3.6)

- In
  [`graph_delta_area()`](https://alinetalhouk.github.io/diceR/reference/graphs.md),
  use base solution instead of
  [`flux::auc()`](https://rdrr.io/pkg/flux/man/auc.html). `flux` also
  depends on `caTools`

- In
  [`prepare_data()`](https://alinetalhouk.github.io/diceR/reference/prepare_data.md),
  use own implementation of `quantable::robustscale()` with all of the
  former function’s defaults. `quantable` also depends on `caTools`

- Specify Bioconductor installation on Travis and AppVeyor since `NMF`
  now Imports `Biobase`

### Minor improvements and bug fixes

- Remove `suppressWarnings(RNGversion("3.5.0"))` after updating R
  version

- Run [`LCA()`](https://alinetalhouk.github.io/diceR/reference/LCA.md)
  unit test on imputed clustering object

- Remove internal validity measures with any `Inf` entries for
  `consensus_reweigh()`

- Use a cleaner, more robust method of removing `Rplots.pdf` after
  running `test-graphs.R`

- Ensure column binding with
  [`purrr::map_dfc()`](https://purrr.tidyverse.org/reference/map_dfr.html)
  in `consensus_rank()`

- Replaced
  [`dplyr::bind_cols()`](https://dplyr.tidyverse.org/reference/bind_cols.html)
  with
  [`purrr::flatten_dfc()`](https://purrr.tidyverse.org/reference/flatten.html)
  to suppress warning “Outer names are only allowed for unnamed scalar”
  in `get_cdf()`

- update roxygen and docs

## diceR 0.6.0

CRAN release: 2019-07-25

- Remove deprecated `dplyr` functions and use `.data` pronoun

- k-means clustering should not support distance matrices as input
  ([@jerryji1993](https://github.com/jerryji1993),
  [\#139](https://github.com/AlineTalhouk/diceR/issues/139))

- Add LCA as a consensus function
  ([@philstraforelli](https://github.com/philstraforelli),
  [\#137](https://github.com/AlineTalhouk/diceR/issues/137))

## diceR 0.5.2

CRAN release: 2019-03-08

- Fix `length > 1 in coercion to logical` error in
  [`consensus_evaluate()`](https://alinetalhouk.github.io/diceR/reference/consensus_evaluate.md)
  due to comparisons using `||` operator

- Add `suppressWarnings(RNGversion("3.5.0"))` before call to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) in examples, tests,
  and vignette to use old RNG sampling

- Use `.covrignore` to exclude `zzz.R` from being considered in code
  coverage

- Use `dplyr` version \>= 0.7.5 to ensure
  [`bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html)
  works

- Fixed bug where scaled matrix using the “robust” method in
  [`prepare_data()`](https://alinetalhouk.github.io/diceR/reference/prepare_data.md)
  was nested in the `data` element
  ([@AlineTalhouk](https://github.com/AlineTalhouk),
  [\#134](https://github.com/AlineTalhouk/diceR/issues/134))

## diceR 0.5.1

CRAN release: 2018-06-11

- Add parameter `hc.method` in `dice` and `consensus_cluster` to pass to
  `method` parameter in
  [`stats::hclust`](https://rdrr.io/r/stats/hclust.html)
  ([@JakeNel28](https://github.com/JakeNel28),
  [\#130](https://github.com/AlineTalhouk/diceR/issues/130))

- Remove dependencies on `largeVis`: package will be archived

## diceR 0.5.0

CRAN release: 2018-05-05

- Revert back to using `NMF` since `NNLM` has been archived and `NMF` is
  back in active maintenance.

- Choose fuzzifier m in `cmeans` using Equation 5 from
  <https://academic.oup.com/bioinformatics/article/26/22/2841/227572>
  (thanks [@Asduveneck](https://github.com/Asduveneck))

## diceR 0.4.0

CRAN release: 2018-02-22

- Replace all code that depended on `NMF` with `NNLM` and `pheatmap`:
  CRAN notified that `NMF` will be archived because of inactive
  maintenance

- Update `.yml` files default templates

## diceR 0.3.2

CRAN release: 2018-01-14

- Fix bug in
  [`consensus_cluster()`](https://alinetalhouk.github.io/diceR/reference/consensus_cluster.md)
  when custom algorithms were excluded from output (thanks
  [@phiala](https://github.com/phiala))

- Use markdown language for documentation

- Various performance improvements and code simplifications

## diceR 0.3.1

CRAN release: 2017-12-12

- Suppress success/fail message printout and fix input data to be matrix
  for block clustering

- Fix bug in `algii_heatmap()` when `k.method = "all"` in
  [`dice()`](https://alinetalhouk.github.io/diceR/reference/dice.md)

- Fix bug in calculating internal indices when data has categorical
  variables (thanks Kurt Salmela)

## diceR 0.3.0

CRAN release: 2017-11-29

- Updated object output names in
  [`consensus_evaluate()`](https://alinetalhouk.github.io/diceR/reference/consensus_evaluate.md)

- Fix unit test in `test-dice.R` for R-devel

- Add internal function: ranked algorithms vs internal validity indices
  heatmap graph

- Fix bugs in
  [`graph_cdf()`](https://alinetalhouk.github.io/diceR/reference/graphs.md),
  [`graph_tracking()`](https://alinetalhouk.github.io/diceR/reference/graphs.md)
  when only one k selected

- Progress messages in
  [`dice()`](https://alinetalhouk.github.io/diceR/reference/dice.md)

- Fix bug in
  [`consensus_evaluate()`](https://alinetalhouk.github.io/diceR/reference/consensus_evaluate.md)
  when algorithm has `NA` for all PAC values

## diceR 0.2.0

CRAN release: 2017-09-29

- New dimension reduction methods: t-SNE, largeVis
  ([@dustin21](https://github.com/dustin21))

- Better annotated progress bar using `progress` package

- Speed up the operation that transforms a matrix to become “NMF-ready”

- Simplify saving mechanism in
  [`consensus_cluster()`](https://alinetalhouk.github.io/diceR/reference/consensus_cluster.md)
  such that only `file.name` needs to be specified, and the `save`
  parameter has been removed

- New algorithms: SOM, Fuzzy C-Means, DBSCAN
  ([@dustin21](https://github.com/dustin21),
  [\#118](https://github.com/AlineTalhouk/diceR/issues/118))

- Added significance testing section to vignette

- Fixed direction of optimization: compactness and connectivity should
  be minimized

## diceR 0.1.0

CRAN release: 2017-06-21

- New submission to CRAN accepted on June 21, 2017
