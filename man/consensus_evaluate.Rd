% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/consensus_evaluate.R
\name{consensus_evaluate}
\alias{consensus_evaluate}
\title{Evaluate, trim, and reweigh algorithms}
\usage{
consensus_evaluate(
  data,
  ...,
  cons.cl = NULL,
  ref.cl = NULL,
  k.method = NULL,
  plot = FALSE,
  trim = FALSE,
  reweigh = FALSE,
  n = 5,
  lower = 0,
  upper = 1
)
}
\arguments{
\item{data}{data matrix with rows as samples and columns as variables}

\item{...}{any number of objects outputted from \code{\link[=consensus_cluster]{consensus_cluster()}}}

\item{cons.cl}{matrix of cluster assignments from consensus functions such as
\code{kmodes} and \code{majority_voting}}

\item{ref.cl}{reference class}

\item{k.method}{determines the method to choose k when no reference class is
given. If \code{ref.cl} is not \code{NULL}, this is the number of distinct classes in
the reference; otherwise the chosen k is determined by the one giving the
largest mean PAC across algorithms. Alternatively, specifying an integer
will override the best chosen k, and specifying "all" will produce
consensus results for all k values ("all" is implicitly used when there is
only one k).}

\item{plot}{logical; if \code{TRUE}, \code{graph_all} is called}

\item{trim}{logical; if \code{TRUE}, algorithms that score low on internal indices
will be trimmed out}

\item{reweigh}{logical; if \code{TRUE}, after trimming out poor performing
algorithms, each algorithm is reweighed depending on its internal indices.}

\item{n}{an integer specifying the top \code{n} algorithms to keep after trimming
off the poor performing ones using Rank Aggregation. If the total number of
algorithms is less than \code{n} no trimming is done.}

\item{lower}{the lower bound that determines what is ambiguous}

\item{upper}{the upper bound that determines what is ambiguous}
}
\value{
\code{consensus_evaluate} returns a list with the following elements
\itemize{
\item \code{k}: if \code{ref.cl} is not \code{NULL}, this is the number of distinct classes
in the reference; otherwise the chosen \code{k} is determined by the one giving
the largest mean PAC across algorithms
\item \code{pac}: a data frame showing the PAC for each combination of algorithm
and cluster size
\item \code{ii}: a list of data frames for all k showing internal evaluation
indices
\item \code{ei}: a data frame showing external evaluation indices for \code{k}
\item \code{trim.obj}: A list with 4 elements
\itemize{
\item \code{alg.keep}: algorithms kept
\item \code{alg.remove}: algorithms removed
\item \code{rank.matrix}: a matrix of ranked algorithms for every internal
evaluation index
\item \code{top.list}: final order of ranked algorithms
\item \code{E.new}: A new version of a \code{consensus_cluster} data object
}
}
}
\description{
Evaluates algorithms on internal/external validation indices. Poor performing
algorithms can be trimmed from the ensemble. The remaining algorithms can be
given weights before use in consensus functions.
}
\details{
This function always returns internal indices. If \code{ref.cl} is not \code{NULL},
external indices are additionally shown. Relevant graphical displays are also
outputted. Algorithms are ranked across internal indices using Rank
Aggregation. Only the top \code{n} algorithms are kept, the rest are trimmed.
}
\examples{
\dontshow{if (rlang::is_installed("apcluster")) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
# Consensus clustering for multiple algorithms
set.seed(911)
x <- matrix(rnorm(500), ncol = 10)
CC <- consensus_cluster(x, nk = 3:4, reps = 10, algorithms = c("ap", "km"),
progress = FALSE)

# Evaluate algorithms on internal/external indices and trim algorithms:
# remove those ranking low on internal indices
set.seed(1)
ref.cl <- sample(1:4, 50, replace = TRUE)
z <- consensus_evaluate(x, CC, ref.cl = ref.cl, n = 1, trim = TRUE)
str(z, max.level = 2)
\dontshow{\}) # examplesIf}
}
