#' Title function to update U, I, and J for function linkage
#'
#' @param i : see i in linkage.R
#' @param j : see j in linkage.R
#' @param m : see m in linkage.R
#'
#' @return a list of U, I, and J in linkage.R
#' @export
#'
#' @examples
UIJ<-function(i,j,m){
  I1<-NULL
  I2<-NULL
  I3<-NULL
  I<-NULL
  J<-NULL
  if(i>=2){
    I1<-1:(i-1)
  }
  if((j-1)>=(i+1)){
    I2<-(i+1):(j-1)
  }
  if(m>=(j+1)){
    I3<-(j+1):m
  }
  if(!is.null(I1)){
    I<-c(I,I1*(m-(I1+1)/2)-m+i)
    J<-c(J,I1*(m-(I1+1)/2)-m+j)
  }
  if(!is.null(I2)){
    I<-c(I,i*(m-(i+1)/2)-m+I2)
    J<-c(J,I2*(m-(I2+1)/2)-m+j)
  }
  if(!is.null(I3)){
    I<-c(I,i*(m-(i+1)/2)-m+I3)
    J<-c(J,j*(m-(j+1)/2)-m+I3)
  }
  U<-c(I1,I2,I3)
  return(list(U=U,I=I,J=J))
}
