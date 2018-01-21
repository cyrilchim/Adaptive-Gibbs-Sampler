#ifndef update_matrix_hpp
#define update_matrix_hpp
// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>
#include "includes/Adaptive_Gibbs.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;


// UPDATE COVARIANCE INVERSE COVARIANCE 'Q' AND THE ITS CHOLESKY DECOMPOSITION 'L_1'
// also we compute a set of matrices S_cond(i), where for each i S_cond(i) = Q^{-1}_ii
// is an inverse of a i-th diagonal block of Q. If the target distribution is normal,
// S_cond is a conditional covariance matrix of the block i

void update_matrix(mat S, mat* Q, mat* L_1, field<mat>* S_cond, param* alg_param)
{
    
    uvec blocking_structure = alg_param->blocking_structure;
    
    int i;
    if(par>1)
    {
        // compute Q^{-1} only if `estimate_spectral_gap` is enabled or some of the blocks are high dimensional
        if(alg_param->estimate_spectral_gap == 1 || par!=dim)
        {
            //reestimation of inverse of the covariance
            if(alg_param->perturb_covariance==1)
            {
                inv_sympd(*Q, S + alg_param->perturb_coeff * eye<mat>(dim,dim));
            }
            else
            {
                inv_sympd(*Q, S);
                if (Q->is_empty()){inv_sympd(*Q, S + lowerb*eye<mat>(dim,dim));}
            }
            
            for(i = 0; i<par; i++)
            {
                
                L_1->submat(blocking_structure[i],blocking_structure[i],blocking_structure[i+1]-1,blocking_structure[i+1]-1)
                = chol(Q->submat(blocking_structure[i],blocking_structure[i],blocking_structure[i+1]-1,blocking_structure[i+1]-1), "lower"); //update diagonal block matrix D^{-1} from the algorithm
                
                if(blocking_structure[i + 1] - blocking_structure[i] > 1)
                {
                    (*S_cond)(i) = inv(L_1->submat(blocking_structure[i], blocking_structure[i], blocking_structure[i+1]-1, blocking_structure[i+1]-1)).t() ;
                }
                
            }
        }

    }
    else
    {
        if(dim > 1)
        {
          chol((*S_cond)(0), S, "lower");
          if ((*S_cond)(0).is_empty()){(*S_cond)(0) = chol(S + lowerb*eye<mat>(dim,dim), "lower");}
        }
    }
    
    
}


#endif /* update_matrix_hpp */