#' Proportion of Ambiguous Clustering
#'
#' Given a consensus matrix, returns the proportion of ambiguous clusters (PAC).
#' This is a robust way to assess clustering performance.
#'
#' Since a consensus matrix is symmetric, we only look at its lower (or upper)
#' triangular matrix. The proportion of entries strictly between \code{lower}
#' and \code{upper} is the PAC. In a perfect clustering, the consensus matrix
#' would consist of only 0s and 1s, and the PAC assessed on the (0, 1) interval
#' would have a perfect score of 0. Using a (0.1, 0.9) interval for defining
#' ambiguity is common as well.
#'
#' The PAC is not, strictly speaking, an internal validity index. Originally
#' used to choose the optimal number of clusters, here we use it to assess
#' cluster stability. However, PAC is still agnostic any gold standard
#' clustering result so we use it like an internal validity index.
#'
#' @param cm consensus matrix. Should be symmetric and values between 0
#'   and 1.
#' @param lower the lower bound that determines what is ambiguous
#' @param upper the upper bound that determines what is ambiguous
#' @return the PAC is a score used in clustering performance. The lower it is
#'   the better, because we want minimal ambiguity amongst the consensus.
#' @author Derek Chiu
#' @references Senbabaoglu, Y., Michailidis, G., & Li, J. Z. (2014). Critical
#'   limitations of consensus clustering in class discovery. Scientific reports,
#'   4.
#' @export
#' @examples
#' set.seed(1)
#' x <- replicate(100, rbinom(100, 4, 0.2))
#' y <- consensus_matrix(x)
#' PAC(y, lower = 0.05, upper = 0.95)
PAC <- function(cm, lower = 0, upper = 1) {
  if (all(cm == 0)) {
    pac <- NA
  } else {
    pac <- cm %>%
      magrittr::extract(lower.tri(.)) %>%
      magrittr::extract(. > lower & . < upper) %>%
      length() %>%
      magrittr::divide_by(., length(cm[lower.tri(cm)]))
  }
  pac
}
