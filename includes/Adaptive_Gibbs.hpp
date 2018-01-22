#ifndef Adaptive_Gibbs_hpp
#define Adaptive_Gibbs_hpp
#define ARMA_DONT_PRINT_ERRORS
// [[Rcpp::depends(RcppArmadillo)]]
#include <cstdlib>
#include <RcppArmadillo.h>

#include <fstream>

using namespace Rcpp;
using namespace arma;
using namespace std;

// GLOBAL VARIABLES

int dim; // dimensionality of the target distribution 
int par; // number of parameters/blocks to update by the Gibbs/MwG algorithms.
         // If no blocking is used, par = d
int burn_in; // number of burn-in samples
double lowerb;// lower bound \epsilon from the ARSG algorithm. By defualt it is set to be 0.1/par^2
double am; // sequence a_m from the ARSG algorithm
double m_iter; // index m of the sequence a_m

uvec blocking_structure; // par+1 dimensional vector that allows easy access to the elements of blocks for the blocking parameter of AMCMC(...) function.  Defined as blocking_structure[0] = 0, blocking_structure[i] - blocking_structure[i-1] = blocking[i-1], i=0 ,.., par

vec current_location; // the current location of the chain

string working_directory; // write/read directory. You can set the directory using set_working_direcoty(directory) function described below


int err_type; // error variable to be used in the constructor of inherited classes of distribution_class. See gaussian_target.hpp for an example of use.  Non-zero value halts the program 
string err_name; // error reference name



typedef Rcpp::Nullable<Rcpp::Function> nullable_func;

// variable bearing information about algorithm aprameters
struct param{
    param() : save(0), burn_in(1), logdensity(0), full_cond(0), gibbs_sampling(0),
              estimate_spectral_gap(1), beta(0), frequency_proposal_update(10) {}
    double *p; //probability weigths for the RSGS algorithm
    bool adapt_proposals;//enable/disable adaptation of proposal varianes in the Metropolis part of the algorithm
    bool adapt_weights;//enable/disable adaptation of sampling weights
    nullable_func R_density;//R-defined density
    bool c_density_flag; //use a density defined in c++
    bool save;//enable/disable saving of the output of the algorithm
    bool burn_in; //burn in mode: 0 - disabled; 1 - enabled;
    int batch_length; //number of samples between adaptation of the sampling weigths
    int thin; //thinning
    bool logdensity; //0 - use density; 1 - use logdensity
    bool full_cond; // 0 - do not use full conditionals to compute acceptance ratio; 1 - use full conditionals
    bool gibbs_sampling; //0 - Metropolis within Gibbs; 1 - Gibbs Samoling
    bool estimate_spectral_gap; // 0 - do not estimate spectral gap; 1 - estimate spectral gap
    uvec gibbs_step; // vector of 0-1 indicating whether block i is updated using Gibbs Sammpler or Metropolis-Hastings normal proposal
    bool ignore_gibbs_sampling; //0 - ignore gibbs_sampling parameter
    uvec blocking_structure; // vector describing blocking structure for the algorithm. The same as blocking_structure
    double beta; //parameter beta>=0 from the AIRMCMC paper that corresponds to the diminishing adaptation rate
    unsigned int frequency_proposal_update;//relative frequency of the proposal variances adaptations. By default, frequency_ratio  = 10, meaning the proposals are updated 10 times more frequently than the weights.
    double current_log_density;//value of  log-density at the current state of the chain
    double proposal_log_density;//value of  log-density at the Metropolis proposal
    bool reset_log_density; // whether to recompute current_log_density
    bool perturb_covariance; // perturb covarince estimation by adding S = S + perturb_coeff * I)
    double perturb_coeff; // perturbation coefficient used by `perturb_covariance`
    double stabilizing_weight; // used to generate stabilise the proposal:
           // proposal = stabilizing_weight * N(current_state, stabilizing_coeff^2) + (1 - stabilising_weight) * adaptive_proposal
    double stabilizing_coeff; // used by for generating a proposal, see `stabilizing_weight`
    
    void print_param ();//print values of a param variable. Should be used for developing/debugging purposes only
};

//CLASS FOR THE TARGET DISTRIBUTION
class distribution_class{
public:
  virtual double density(vec theta) = 0;    //target density at the point theta 
  virtual double logdensity(vec theta) = 0; //logarithm of the target density
  
  virtual double full_cond(vec theta, int gp) = 0;  //full conditional of the block gp
  virtual double logfull_cond(vec theta, int gp)  = 0;  //logarithm of the full conditional of the block gp
  
  virtual vec sample_full_cond(vec theta, int gp) = 0;   //sample block gp from its full conditional distribution
  virtual ~distribution_class(){};
};

/*****************************************************************************************/


#endif /* Adaptive_Gibbs_hpp */
