#ifndef adaptive_step_hpp
#define adaptive_step_hpp
// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>
#include "Adaptive_Gibbs.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;


//ADAPTIVE STEP OF THE ARSGS
// perturb - gaussian noise that perturbs the maximum eigenvector
// mu - estimated mean of the target distribution
// tr_aux - auxiliary variable that records samples that are used by adaptive_step in order to update the covariance S. This variable is deleted after the call of adaptive_step
// start_old - last index of a sample that was already used in updating the covariance S
// start_new - new value for start_old after the samples of tr_aux are used to update S
// dir - gradient direction from ARSG algorithm
// w, p - probability weights from ARSG algorithm
// inv_sp_gap, inv_sp_gap_trace - estimated spectral gap and its trace
// tr_proposals - trace of the scale parameters of Metropolis part of the algorithm

void adaptive_step(vec *perturb, mat *S, vec* mu, int* start_old, int* start_new, vector<vector<double> >* tr_aux,  mat* Q, mat* L_1, field<mat>* S_cond, vec* scale, vec* dir, vec* z, vec* w, vec* p, double* inv_sp_gap, vector<double> *inv_sp_gap_trace, vector<vector<double> > *tr_proposals, param *alg_param)
{
    
    double w_array[par+1]; //array version of the weigths w
    double p_array[par]; //array version of the weigths p
    int i, j;
    vec temp_mu(dim);
    vec theta(dim);
    
    
    //COVARIANCE MATRIX REESTIMATION
    if(alg_param->adapt_weights==1 || alg_param->estimate_spectral_gap==1 || alg_param->adapt_proposals==1)
    {
        j = *start_old + 1;
        while(j <= *start_new)
        {
            for(i = 0; i<dim; i++)
            {
                theta(i) = (*tr_aux)[i].back();//get last element of the auxiliary vector (*tr_aux)[i]
                (*tr_aux)[i].pop_back(); //erase last element
            }
            
            
            if((alg_param->adapt_weights==1||alg_param->estimate_spectral_gap==1||alg_param->adapt_proposals==1)&&(alg_param->burn_in==0))
            {
                //update the estimated covariance matrix
                temp_mu = 1.0 / (j + 1)*(j * (*mu) + theta);
                *S = (j - 1.0) / ( j )*(*S) + (1.0 / (j)*(theta))*(theta).t() + (*mu)*((*mu).t()) - ((j  + 1.0) / (j)*temp_mu)*temp_mu.t();
                *mu = temp_mu;
            }
            j = j + 1;
        }
        
        update_matrix(*S, Q, L_1, S_cond, alg_param);//Step 2 of the ARSGS
        
    }
    
    *start_old = *start_new;
    
    
    
    if(alg_param->estimate_spectral_gap == 1||alg_param->adapt_weights == 1)
    {
        for(i=0;i<par;i++) {
            p_array[i] = (*p)(i);  //transform vec to array
            w_array[i] = (*w)(i) ; //transform vec to array
        }
        w_array[par] = (*w)(par) ; //transform vec to array
        
        
        *dir = grad_direction(perturb, z, S, L_1, w_array, inv_sp_gap, alg_param);
        if (alg_param->adapt_weights == 1) {update_weights(w_array, p_array, *dir);}
        inv_sp_gap_trace->push_back(*inv_sp_gap);
        
        for(i=0;i<par;i++) {
            (*p)(i) = p_array[i]; //transform array back to vec
            (*w)(i) = w_array[i]; //transform array back to vec
        }
        (*w)(par) = w_array[par];
        
        
    }
    
    
    
    //record the trace of the scale parameters
    if(alg_param->adapt_proposals == 1)
    {
        for(int i=0;i<par;i++)
        {
            (*tr_proposals)[i].push_back( (*scale)(i) );
        }
    }
    
    
}


#endif /* adaptive_step_hpp */