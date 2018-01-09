#ifndef poisson_hierarchical_model_hpp
#define poisson_hierarchical_model_hpp
#define ARMA_DONT_PRINT_ERRORS
// [[Rcpp::depends(RcppArmadillo)]]
#include "includes/Adaptive_Gibbs.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;



// FEEL FREE TO USE BUT DO NOT CHANGE THE GLOBAL VARIABLES IN THIS FILE (THE GLOBAL VARIABLES ARE LISTED IN Adaptive_Gibbs.hpp)

class poisson_hierarchical_model: public distribution_class{
public:
    mat X; // design matrix
    vec Y; // data vector
    vec aux; // auxiliary vector for fast density computation
    int n_data; // number of data points
    int last_ind; // last updated block
    int call_num; // number of density estimations
    vec prior_mean, prior_var; // mean and variances of the prior distribution
    vec posterior_part; //part of the mean of the posterior
    // add some additional functions/parameters if necessary
    poisson_hierarchical_model();//specify constructor  if necessary
    ~poisson_hierarchical_model(); //specify destructor if necessary
    double density(vec theta);    //target density at the point theta 
    double logdensity(vec theta); //logarithm of the target density
    
    double full_cond(vec theta, int gp);  //full conditional of the block gp
    double logfull_cond(vec theta, int gp);  //logarithm of the full conditional of the block gp
    vec sample_full_cond(vec theta, int gp);  //sample block gp from its full conditional distributiont
};


poisson_hierarchical_model::poisson_hierarchical_model()
{
  call_num = 0;
  last_ind = 0;
  
  prior_mean.set_size(dim);
  prior_mean.ones();
  prior_mean = (-1) * prior_mean;
  
  prior_var.set_size(dim);
  prior_var.ones();
  
  posterior_part.set_size(dim);
  
  if(X.load(working_directory+"phm_design")==false || Y.load(working_directory+"phm_observations")==false)
  {
    err_name = "error: can't initialise gaussian distribution: file not found";
    err_type = 2;
  }
  else{
    X.load(working_directory+"phm_design");
    Y.load(working_directory+"phm_observations");
    n_data = Y.n_elem;
    
    if(X.n_cols!=dim || X.n_rows!=n_data)
    {
      err_name = "error: can't initialise gaussian distribution: wrong dimensionality or the input";
      err_type = 2;
    }
  }
  
  for (int i = 0; i<dim; i++) {
    posterior_part(i) = prior_mean(i) / prior_var(i);
    for (int j = 0; j<n_data; j++) {
      posterior_part(i) = posterior_part(i) + X(j, i)*Y(j);
    }
    
  }
  
  // set up the auxiliary vector
  aux.set_size(n_data);
  for(int i=0; i<n_data; i++) //recompute c[i]
  {
    aux[i] = 0;
    for(int j=0; j<dim; j++) { aux[i] = aux[i] + X(i, j)* current_location(j); }
  }
  
  
}

poisson_hierarchical_model::~poisson_hierarchical_model()
{ 
  
} 



double poisson_hierarchical_model::density(vec theta)
{
  return 0;
}

double poisson_hierarchical_model::logdensity(vec theta)
{
  return 0;
}

double poisson_hierarchical_model::full_cond(vec theta, int prob_ind)
{
  return 0;
}


double poisson_hierarchical_model::logfull_cond(vec theta, int prob_ind)
{
  double res;
  

  for(int i=0; i<n_data; i++)
  {
    aux[i] = aux[i] - X(i, prob_ind)*current_location[prob_ind];
  }
  if(call_num>0)
  {
    for(int i=0; i<n_data; i++)
    {
      aux[i] = aux[i] +  X(i, last_ind)*current_location[last_ind];
    }
  }
  
  
  res = -1. / (2 * prior_var(prob_ind)) * theta(prob_ind) * theta(prob_ind) + posterior_part(prob_ind) * theta(prob_ind);
  
  for (int i=0; i<n_data; i++) 
  {
    res = res - exp(aux[i] + X(i, prob_ind)*theta(prob_ind));     
  }
  
  call_num++;
  last_ind = prob_ind;
  return res;
}


//should return a vector of the same dimension as the block gp
vec poisson_hierarchical_model::sample_full_cond(vec theta, int prob_ind)
{
  vec res;
  return res;
}


//Precompiled C++ functions for use in R 

// [[Rcpp::export]]
void phm_example_1_data(int N_data, int N_param) // create matrix X for large correlations
{
  mat X(N_data, N_param); // design matrix
  X.zeros();
  
  int i, j;
  
  int k = int(N_data/N_param);
  
  for(i=0;i<2*k;i++){
    X(i,0) = 1;
    X(i,1) = 1;
  }
  
  for(i=2*k;i<5*k;i++){ 
    X(i,2) = 1;
    X(i,3) = 1;
    X(i,4) = 1;
  }
  
  
  
  for (j=5; j<N_param; j++) {
    for (i=j*k; i<(j+1)*k; i++) {
      X(i,j) = 1;
    }
  }
  
  
  
  for (i = 0; i<N_data; i++) {         // perturb the design matrix
    for (j = 0; j<N_param; j++) {
      X(i,j) = X(i,j) + 0.1 * rbeta(1, 0.1, 0.1)[0];
    }
  }
  
  X.save(working_directory+"phm_design");
  
  vec beta_true(N_param);// true `beta` parameter
  beta_true.ones(); // set `beta` to (1,..,1)
  
  vec lambda_true(N_data); // true value for lambda
  vec Y(N_data); // observations
  
  for (j = 0; j<N_data; j++) 
  {
    lambda_true(j) = 1;
    for (i = 0; i<N_param; i++)
    {
      lambda_true(j) = lambda_true(j)* exp(X(j, i)*beta_true(i));
    }
  }
  for (j = 0; j<N_data; j++) {
    Y(j) = rpois(1, lambda_true(j))[0]; //generate data Y with seed r0
  }
  
  Y.save(working_directory+"phm_observations");
}

// [[Rcpp::export]]
NumericMatrix phm_get_design_mat()
{
  mat X;
  X.load(working_directory+"phm_design");
  return wrap(X);
}

// [[Rcpp::export]]
NumericVector phm_get_observations()
{
  vec Y;
  Y.load(working_directory+"phm_observations");
  return wrap(Y);
}
#endif /* poisson_hierarchical_model_hpp */
