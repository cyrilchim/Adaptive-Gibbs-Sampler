#ifndef multinomial_hpp
#define multinomial_hpp
// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;


//sample from a categorical distribution with probability weights p
int multinomial(double *p) {
    int res=0;
    double u = runif(1)[0];
    double cumsum = p[res];
    
    
    while(u>cumsum && (res < par))
    {
        res++;
        if(res < par - 1){cumsum += p[res];}
        else {cumsum = 1.;} // the condition prevents crashing in case \sum p_i \neq 1
    }
    
    return res;
}


#endif /* multinomial_hpp */