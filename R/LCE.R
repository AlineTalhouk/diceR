#' Generate CTS, SRS, and ASRS similarity matrices.
#' 
#' 
#' @param E is a 2D matrix of clustering results. Complete and relabelled cases needed.
#' @param dcCTS decay constant for CTS matrix
#' @param dcSRS decay constant for SRS matrix
#' @param dcASRS decay constant for ASRS matrix
#' @param R number of repetitions for SRS matrix
#' @author Johnson Liu
#' @return a list containing the CTS, SRS, and ASRS matrix
#' @importFrom stats complete.cases
#' @export
#' 
#' @examples
#' set.seed(1)
#' E <- matrix(rep(sample(1:4, 1000, replace = TRUE)), nrow = 100, byrow = FALSE)
#' link_clust(E = E, dcCTS = 0.8, dcSRS = 0.8, dcASRS = 0.8, R = 10,
#'            is.relabelled = FALSE)
#' data(hgsc)
#' dat <- t(hgsc[, -1])
#' x <- ConClust(dat[1:100, 1:50], k = 4, reps = 4,
#'               method = c("nmfEucl", "hcAEucl", "hcDianaEucl"), save = FALSE)
#' y <- link_clust(E = x, dcCTS = 0.8, dcSRS = 0.8, dcASRS = 0.8, R = 10)
LCE <- function(E, dcCTS = 0.8, dcSRS = 0.8, dcASRS = 0.8, R = 10,
                similarity_matrices="cts",
                is.relabelled = TRUE) {
  assertthat::assert_that(is.array(E) || is.matrix(E), dcCTS >= 0 && dcCTS <= 1,
                          dcASRS >= 0 && dcASRS <= 1, dcSRS >= 0 && dcSRS <= 1,
                          is_pos_int(R), sum(!similarity_matrices%in%c("cts","asrs","srs"))==0)
  CTS=cts(E=E,dc=dcCTS)
  ASRS<-NULL
  SRS<-NULL
  if("asrs"%in%similarity_matrices){
    ASRS<-asrs(E=E,dc=dcASRS)
  }
  if("srs"%in%similarity_matrices){
    SRS<-srs(E=E,dc=dcSRS,R=R)
  }
  return(list(CTS=CTS,SRS=SRS,ASRS=ASRS))
}
