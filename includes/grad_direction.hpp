#ifndef grad_direction_hpp
#define grad_direction_hpp
// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>
#include "Adaptive_Gibbs.hpp"

using namespace Rcpp;
using namespace arma;


//APROXIMATE GRADIENT DIRECTION

vec grad_direction(vec *perturb, vec* z, mat* S, mat *L_1, double* w, double* inv_sp_gap, param *alg_param)
{
    
    vec z1(dim+1), dir(par);
    int i;
    double s=0;
    z1.zeros();
    
    uvec blocking_structure = alg_param->blocking_structure;
    
    mat L(dim, dim);
    L.zeros();
    
    for(i=0;i<par; i++)
    {
        L.submat(blocking_structure[i],blocking_structure[i],blocking_structure[i+1]-1,blocking_structure[i+1]-1) = 1./sqrt(w[i])*L_1->submat(blocking_structure[i],blocking_structure[i],blocking_structure[i+1]-1,blocking_structure[i+1]-1);
    }
    
    
    //for(l=0; l<npi; l++){
    z1.zeros();
    
    if(alg_param->perturb_covariance==1)
    {
        z1.subvec(0,dim-1) = L.t() * ((*S  + 1. / (par * par * par) * eye<mat>(dim,dim)) * (L * (*z).subvec(0,dim-1) ));
    }
    else
    {
        z1.subvec(0,dim-1) = L.t() * ((*S) * (L * (*z).subvec(0,dim-1) ));
    }
    
    z1(dim) = 1./w[par] * (*z)(dim);
    
    //  if (l == npi - 1) {
    for(i=0; i<par; i++) {s = s + w[i];}
    *inv_sp_gap =  s * norm(z1, 2) / norm(*z, 2);
    // }
    *z = z1;
    *z = normalise(*z); // normalization
    //  }
    
    *z = *z + sqrt(am) * normalise(*perturb);//perturbation to improve stability
    
    //"descend" direction
    *z = normalise(*z); // normalisation
    
    
    for(i=0;i<par;i++)
    {
        dir(i) =((1.0 / w[i])*dot((*z).subvec(blocking_structure[i],blocking_structure[i+1]-1),(*z).subvec(blocking_structure[i],blocking_structure[i+1]-1)) - (1.0 / w[par])*(*z)(dim)*(*z)(dim));
    }
    
    dir = normalise(dir,1); 
    
    return(dir);
    
}


#endif /* grad_direction_hpp */