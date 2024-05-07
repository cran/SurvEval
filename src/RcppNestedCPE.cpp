#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double reducedCPE(
    NumericMatrix bMatrix, 
    NumericMatrix gMatrix, 
    NumericMatrix bKernel, 
    NumericMatrix SurvMatrix
){
  
  int na =  bMatrix.nrow(); 
  
  int i, j, k, l;
  
  double numcounter;
  double denomcounter;
  double outsidesum=0;
  
  for(i=0; i< (na-1); i ++){
    for(j = (i+1); j < na; j ++){
      
      numcounter = 0;
      denomcounter = 0;
      for(k=0; k< na; k ++){
        for(l = 0; l < na; l ++){
          
          if(l != k){
            numcounter +=      (1/( 1+exp(bMatrix(i,j)+gMatrix(k,l))))*(1-SurvMatrix(i,k)*SurvMatrix(j,l))*bKernel(i,k)*bKernel(j,l);
            denomcounter +=      bKernel(i,k)*bKernel(j,l);
          }
        }                  
      }
      
      outsidesum += (numcounter/denomcounter); 
    }
  }
  
  return(outsidesum);
}
