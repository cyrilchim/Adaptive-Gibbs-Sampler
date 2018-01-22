#ifndef sample_batch_hpp
#define sample_batch_hpp
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "Adaptive_Gibbs.hpp"

using namespace Rcpp;
using namespace arma;

// acceptance probability in Metropolis algorithm computed from density ratio
double acceptance_density(vec prop, vec theta, int prob_ind,
                          distribution_class* c_density, param* alg_param);

// acceptance probability in Metropolis algorithm computed from log-density ratio
double acceptance_logdensity(vec prop, vec theta, int prob_ind, distribution_class* c_density,
                             param* alg_param);



//RSGS STEP (Step 1 from ARSGS paper)
// here theta is a starting location for a batch, scaling - vector of scalings in Metropolis-within-Gibbs algorithm
// tr - trace of the chain
// tr_aux - auxiliary variable that records samples that are used by adaptive_step in order to update the covariance S
// start - index of the variable from which the batch is sampled
// count - batch number
void sample_batch(vec* theta, vec* scaling, mat *L_1, field<mat>* S_cond, unsigned int* count,
                  vector<vector<double> > *tr, vector<vector<double> > *tr_aux, int *start,
                  distribution_class* c_density, param *alg_param);
#endif /* sample_batch_hpp */