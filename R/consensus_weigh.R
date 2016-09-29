#' @param data.index matrix of indices for each algorthim outputted from
#'   \code{consensus_compare}.
#' @param top how many top algorithms to include for weighting. Defaults to 5.
#' @return \code{consensus_weigh} returns a matrix of weighted algorithms.
#' @rdname consensus_combine
#' @importFrom utils head
#' @export
consensus_weigh <- function(data.index, top = 5) {
  Weight.PAC <- Weight.CHI <- CHI <- NULL
  wts <- data.index %>%
    arrange(desc(CHI)) %>%
    head(top) %>%
    mutate(Weight.PAC = (1 - PAC) / sum((1 - PAC)),
           Weight.CHI = CHI / sum(CHI),
           Weight = (Weight.PAC + Weight.CHI) / 2)
  return(wts)
}
