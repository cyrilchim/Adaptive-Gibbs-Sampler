#ifndef projection_hpp
#define projection_hpp
// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>
#include "Adaptive_Gibbs.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;


//PROJECTION ALGORITHM FROM ARSGS PAPER

void projection(double *p, int k)
{
    
    vec proj(k), u(k-1), proj_1(k-1);
    int i, rho;
    double s=0, s_1, s_2, lambda=1;
    
    for(i=0; i<k-1; i++)
    {
        proj(i) = p[i];
        if(proj(i)<lowerb) {proj(i) = lowerb; }
        s = s + proj(i);
    }
    
    proj(k-1) = 1 - s;
    
    
    if(proj(k-1)>lowerb)
    {
        for(i=0; i<k; i++) p[i] = proj(i);
    }
    
    else
    {
        p[k-1] = lowerb;
        proj_1 =1. / (1 - k * lowerb) * (proj.subvec(0, k-2) - lowerb);
        u = sort(proj_1,"descend");
        
        rho = 0; s_2 = u(0); s_1 = u(0) + 1 - s_2; lambda = s_1;
        for (i = 1; i < k - 1; i++) {
            s_2 = s_2 + u(i); s_1 = u(i) + 1. / (i + 1) * (1 - s_2);
            if(lambda > 0){ lambda = s_1; rho = i;}
        }
        
        lambda = lambda - u(rho);
        for (i = 0; i < k-1; i++) {
            if(proj_1(i) + lambda >0) {p[i] = lowerb + (1-k*lowerb)*(proj_1(i)+lambda);}
            else {p[i] = lowerb;}
        }
        
    }
    
    
}



#endif /* projection_hpp */