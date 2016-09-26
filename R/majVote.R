#' Title function to perform major vote on a 3D array with NA entries
#'
#' @param arr : input character array with NAs
#'
#' @return a 3D numeric array with NAs in original array voted
#' @export
#'
#' @examples
majVote<-function(arr){
  arr_out<-array(dim=dim(arr))
  for(k in 1:dim(arr_out)[3]){
    temp<-arr[,,k]
    if(sum(!complete.cases(temp))==0){
      arr_out[,,k]<-apply(temp,c(1,2),as.numeric)
    }else{
      for(i in 1:nrow(temp)){
        mode<-names(sort(-table(temp[i,])))[1]
        if(is.na(mode)){
          stop("all entries in a row in array for relabelled imputed array are NAs, you probably need to increase repetitions")
        }
        temp[i,which(!complete.cases(temp[i,]))]<-mode
      }
      arr_out[,,k]<-apply(temp,c(1,2),as.numeric)
    }
  }
  return(arr_out)
}
