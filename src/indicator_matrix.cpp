#include <Rcpp.h>
using namespace Rcpp;

// Indicator matrix
// @noRd
// [[Rcpp::export]]
NumericMatrix indicator_matrix(NumericVector x) {
  int n = x.size();
  NumericMatrix mat(n, n);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      if (!NumericVector::is_na(x[j]) && !NumericVector::is_na(x[i]))
        mat(i, j) = 1;
    }
  }
  return(mat);
}
