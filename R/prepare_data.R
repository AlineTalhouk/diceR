#' Prepare data for consensus clustering
#'
#' Perform feature selection or dimension reduction to remove noise variables.
#'
#' We can apply a basic filtering method of feature selection that removes
#' variables with low signal and (optionally) scales before consensus
#' clustering. Or, we can use t-SNE dimension reduction to transform the data to
#' just two variables. This lower-dimensional embedding allows algorithms such
#' as hierarchical clustering to achieve greater performance.
#'
#' @param data data matrix with rows as samples and columns as variables
#' @param scale logical; should the data be centered and scaled?
#' @param type if we use "conventional" measures (default), then the mean and
#'   standard deviation are used for centering and scaling, respectively. If
#'   "robust" measures are specified, the median and median absolute deviation
#'   (MAD) are used. Alternatively, we can apply "tsne" for dimension reduction.
#' @param min.var minimum variability measure threshold used to filter the
#'   feature space for only highly variable features. Only features with a
#'   minimum variability measure across all samples greater than `min.var` will
#'   be used. If `type = "conventional"`, the standard deviation is the measure
#'   used, and if `type = "robust"`, the MAD is the measure used.
#' @return dataset prepared for usage in `consensus_cluster`
#' @author Derek Chiu
#' @export
#' @examples
#' set.seed(2)
#' x <- replicate(10, rnorm(100))
#' x.prep <- prepare_data(x)
#' dim(x)
#' dim(x.prep)
prepare_data <- function(data, scale = TRUE,
                         type = c("conventional", "robust", "tsne"),
                         min.var = 1) {
  type <- match.arg(type)
  if (type == "tsne") {
    if (!requireNamespace("Rtsne", quietly = TRUE)) {
      stop("Package \"Rtsne\" is needed. Please install it.",
           call. = FALSE)
    } else {
      return(
        data %>%
          as.matrix() %>%
          Rtsne::Rtsne(perplexity = 5, max_iter = 500) %>%
          magrittr::extract2("Y") %>%
          magrittr::set_rownames(rownames(data))
      )
    }
  }
  var.fun <- switch(type, conventional = stats::sd, robust = stats::mad)
  dat <- data %>%
    magrittr::extract(stats::complete.cases(.),
                      apply(., 2, var.fun, na.rm = TRUE) > min.var)
  if (scale) {
    dat <- switch(type,
                  conventional = scale(dat),
                  robust = robust_scale(dat))
  }
  dat
}

#' Same as `quantable::robustscale()` with all default arguments
#' @noRd
robust_scale <- function(data) {
  medians <- apply(data, 2, stats::median, na.rm = TRUE)
  data <- sweep(data, 2, medians, "-")
  mads <- apply(data, 2, stats::mad, na.rm = TRUE)
  mads <- mads / mean(mads)
  data <- sweep(data, 2, mads, "/")
  data
}
