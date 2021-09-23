#include <Rcpp.h>


// [[Rcpp::export]]
Rcpp::NumericVector findEqualGreaterM( Rcpp::NumericVector x, Rcpp::NumericVector values){
  int idx = 0;
  int sizeres = values.size();
  Rcpp::NumericVector res(sizeres);
  
  for(int i = 0; i< sizeres ; i++){
    while ( ( idx < x.size() ) &&  ( x[idx] < values[i] ) ) {
      idx++;
	}
    res[i]=idx+1;
  }
  
  return res;
}





