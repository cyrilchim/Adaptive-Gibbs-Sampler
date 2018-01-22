#ifndef reweight_hpp
#define reweight_hpp
// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>
#include "Adaptive_Gibbs.hpp"

using namespace Rcpp;
using namespace arma;


//reweight probability weights for Metropolis update
void reweight_p(param* alg_param)
{
    double s = 0;
    int i;
    
    for(i = 0; i<par; i++)
    {
        if((alg_param->gibbs_step(i) == 0 && alg_param->ignore_gibbs_sampling == 1) ||
           (alg_param->gibbs_sampling == 0 && alg_param->ignore_gibbs_sampling == 0))
        {
            alg_param->p[i] = (alg_param->blocking_structure(i+1) - alg_param->blocking_structure(i)) * alg_param->p[i];
        }
        s = s + alg_param->p[i];
    }
    
    for(i = 0; i<par; i++)
    {
        alg_param->p[i] = alg_param->p[i] / s;
    }
}

#endif /* reweight_hpp */