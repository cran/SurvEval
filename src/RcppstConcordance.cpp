#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double sstartconcordance(NumericVector Xbetahat, NumericVector Sstart, NumericVector Send){
  int na =  Xbetahat.size(); 
  int i, j;
  
  double output=0; 
  double numcounter = 0;
  double denomcounter = 0;
  
  for(i=0; i< na; i++){
    for(j =0; j < na; j ++){
      if(i != j){
        if(Xbetahat[i] <  Xbetahat[j] ){ 
          numcounter += (1/(1+exp(Xbetahat[i]-Xbetahat[j])))*(Sstart[i]*Sstart[j] - Send[i]*Send[j]); 
        }
        if(Xbetahat[i] !=  Xbetahat[j] ){ 
          denomcounter += (Sstart[i]*Sstart[j] - Send[i]*Send[j]);   
        }                  
      }
    }
  }
  output = (2*numcounter)/denomcounter; 
  return(output);
}

