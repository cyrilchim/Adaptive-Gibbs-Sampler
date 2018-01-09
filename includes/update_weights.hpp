#ifndef update_weights_hpp
#define update_weights_hpp
// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>
#include "Adaptive_Gibbs.hpp"
#include "projection.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;


//UPDATE SAMPLING WEIGHTS (subgradient step of the algorithm together with the projection step)

void update_weights(double* w, double *p, vec dir)
{
    int i, f=1;
    double w1[par+1], s=0;
    for (i = 0; i < par; i++) {
        w1[i] = w[i] + am * dir(i); //subgradient method move
        if (w1[i]<lowerb) { f = 0; }
        s = s + w1[i];
    }
    w1[par] = 1 - s;
    f = f*(w1[par]>=lowerb);
    
    if (f == 1) {
        s = 0;
        for (i = 0; i <= par-1; i++) {
            w[i] = w1[i];
            s = s+w[i];
            
        }
        w[par] = w1[par];
        for (i = 0; i <= par-1; i++) {
            p[i] = 1./s * w[i];
            
        }
        
    }
    else
    {
        for (i = 0; i <= par; i++) w[i] = w1[i];
        
        projection(w, par+1);
        
        s = 0;
        for (i = 0; i <= par-1; i++) {
            
            s = s + w[i];
            
        }
        
        for (i = 0; i <= par-1; i++) {
            p[i] = 1. / s * w[i];
            
        }
        
        
    }
    
    m_iter = m_iter + 1;
    am = log(m_iter) / (m_iter);
}



#endif /* update_weights_hpp */