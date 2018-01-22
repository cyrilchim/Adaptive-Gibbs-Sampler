#ifndef markov_switching_model_hpp
#define markov_switching_model_hpp
#define ARMA_DONT_PRINT_ERRORS
// [[Rcpp::depends(RcppArmadillo)]]
#include "include/Adaptive_Gibbs.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;



//FEEL FREE TO USE BUT DO NOT CHANGE THE GLOBAL VARIABLES IN THIS FILE (THE GLOBAL VARIABLES ARE LISTED IN Adaptive_Gibbs.hpp)

class markov_switching_model: public distribution_class{
public:
    double sw_rate_0, sw_rate_1; // switching rate in regimes 0 and 1, respectively
    double sigma_0, sigma_1; // standard deviations in regimes 0 and 1, respectively
    vec Y; // vector of observation 
    double beta; // observation error

    markov_switching_model(); // specify constructor  if necessary
    ~markov_switching_model(); // specify destructor if necessary
    double density(vec theta); // target density at the point theta 
    double logdensity(vec theta); // logarithm of the target density
    
    double full_cond(vec theta, int prob_ind);  // full conditional of the block gp
    double logfull_cond(vec theta, int prob_ind);  // logarithm of the full conditional of the block gp
    vec sample_full_cond(vec theta, int prob_ind);  // sample block gp from its full conditional distributiont
    
private:
    double sample_regime(vec theta, int prob_ind); // sample conditional regime of the chain
    double sample_state(vec theta, int prob_ind); // sample conditional state of the chain
};


markov_switching_model::markov_switching_model()
{
  sw_rate_0 = 0.001;
  sw_rate_1 = 0.005;
  sigma_0 = sqrt(1.);
  sigma_1 = sqrt(10.);
  beta = sqrt(1.);

  if(Y.load(working_directory+"msm_observations")==0){
    err_name = "error: can't initialise markov switching model: file not found.";
    err_type = 2;
  }
  else{
    Y.load(working_directory+"msm_observations");
    if(Y.n_elem!=dim/2)
    {
      err_name = "error: can't initialise markov switching model: wrong dimensionality of the input";
      err_type = 2;
    }
  }
  
}

markov_switching_model::~markov_switching_model()
{ 
  
} 



double markov_switching_model::density(vec theta)
{
  return 0;
}

double markov_switching_model::logdensity(vec theta)
{
  return 0;
}

double markov_switching_model::full_cond(vec theta, int prob_ind)
{
  return 0;
}


double markov_switching_model::logfull_cond(vec theta, int prob_ind)
{
  return 0;
}


//should return a vector of the same dimension as the block gp
vec markov_switching_model::sample_full_cond(vec theta, int prob_ind)
{
  vec res(1);
  
  if(prob_ind<dim/2)
  {
    res(0) = sample_state(theta, prob_ind);
  }
  else
  {
    res(0) = sample_regime(theta, prob_ind - dim/2);
  }
  
  return res;
}

double markov_switching_model::sample_state(vec theta, int prob_ind)
{
  double res = 0, sd = 1, sd_0 = 1, sd_1 = 1;
  
  vec reg(dim/2), state(dim/2);
  
  for(int i=0; i<dim/2; i++)
  {
    reg(i) = theta(par/2+i);
    state(i) = theta(i);
  }
  

  if(prob_ind==0)
  {
    sd = 1. / (beta * beta);
    if(reg(0)==0) sd = sd + 1. / (sigma_0 * sigma_0);
    else sd = sd + 1. / (sigma_1 * sigma_1);
    
    if(reg(1)==0)
    {
      sd_0 = 1. / (sigma_0 * sigma_0);
      sd = sd + sd_0;
    }
    else {
      sd_0 = 1. / (sigma_1 * sigma_1);
      sd = sd + sd_0;
    }
    sd = 1. / sd;
    res = sd * (state(1) * sd_0 + Y(0) / (beta * beta)) + sqrt(sd) * rnorm(1)[0];
  }
  
  
  
  
  else if (prob_ind==dim/2 - 1)
  {
    sd =1. / (beta * beta);
    if(reg(dim/2 - 1)==0)
    {
      sd_0 =  1. / (sigma_0 * sigma_0);
      sd = sd + sd_0;
    }
    else
    {
      sd_0 =  1. / (sigma_1 * sigma_1);
      sd =  sd + sd_0;
    }
    
    
    sd = 1. / sd;
    res = sd * (state(prob_ind - 1) * sd_0 + Y(prob_ind) / (beta * beta)) + sqrt(sd) * rnorm(1)[0];
    
  }
  

  else
  {
    sd =1. / (beta * beta);
    
    if(reg(prob_ind)==0)
    {
      sd_0 =  1. / (sigma_0 * sigma_0);
      sd = sd + sd_0;
    }
    else
    {
      sd_0 =  1. / (sigma_1 * sigma_1);
      sd =  sd + sd_0;
    }
    
    
    if(reg(prob_ind + 1)==0)
    {
      sd_1 =  1. / (sigma_0 * sigma_0);
      sd = sd + sd_1;
    }
    else
    {
      sd_1 =  1. / (sigma_1 * sigma_1);
      sd =  sd + sd_1;
    }

    sd = 1. / sd;
    res = sd * (state(prob_ind - 1) * sd_0 + state(prob_ind + 1) * sd_1 + Y(prob_ind) / (beta * beta)) + sqrt(sd) * rnorm(1)[0];
    
  }
  
  return res;
}


double markov_switching_model::sample_regime(vec theta, int prob_ind)
  
{
  double res = 0, p_0, p_1, u;
  
  vec reg(dim/2), state(dim/2);

  for(int i=0; i<dim/2; i++)
  {
    reg(i) = theta(par/2+i);
    state(i) = theta(i);
  }

  
  if(prob_ind==0)
  {
    p_0 = 1 / sigma_0 * exp( -(state(prob_ind)) * (state(prob_ind)) / (2 * sigma_0*sigma_0)) * (sw_rate_1 / sw_rate_0) / (1 + sw_rate_1 / sw_rate_0);
    p_1 = 1 / sigma_1 * exp( -(state(prob_ind)) * (state(prob_ind)) / (2 * sigma_1 * sigma_1)) / (1 + sw_rate_1 / sw_rate_0);
    
    
    if(reg(prob_ind + 1)==0)
    {
      p_0 = p_0 * (1 - sw_rate_0);
      p_1 = p_1 * sw_rate_1;
    }
    else
    {
      p_0 = p_0 * sw_rate_0;
      p_1 = p_1 * (1 - sw_rate_1);
    }
    
  }
  
  else if (prob_ind == dim/2-1)
  {
    p_0 = 1 / sigma_0 * exp( -(state(prob_ind) - state(prob_ind - 1)) * (state(prob_ind) - state(prob_ind - 1)) / (2 * sigma_0 * sigma_0));
    p_1 = 1 / sigma_1 * exp( -(state(prob_ind) - state(prob_ind - 1)) * (state(prob_ind) - state(prob_ind - 1)) / (2 * sigma_1 * sigma_1));
    
    if(reg(prob_ind - 1)==0)
    {
      p_0 = p_0 * (1 - sw_rate_0);
      p_1 = p_1 * sw_rate_0;
    }
    else
    {
      p_0 = p_0 * sw_rate_1;
      p_1 = p_1 * (1 - sw_rate_1);
    }
    
    
  }
  

  else
  {
    p_0 = 1 / sigma_0 * exp( -(state(prob_ind) - state(prob_ind - 1)) * (state(prob_ind) - state(prob_ind - 1)) / (2 * sigma_0*sigma_0));
    p_1 = 1 / sigma_1 * exp( -(state(prob_ind) - state(prob_ind - 1)) * (state(prob_ind) - state(prob_ind - 1)) / (2 * sigma_1*sigma_1));
    
    if(reg(prob_ind - 1)==0)
    {
      p_0 = p_0 * (1 - sw_rate_0);
      p_1 = p_1 * sw_rate_0;
    }
    else
    {
      p_0 = p_0 * sw_rate_1;
      p_1 = p_1 * (1 - sw_rate_1);
    }
    
    if(reg(prob_ind + 1)==0)
    {
      p_0 = p_0 * (1 - sw_rate_0);
      p_1 = p_1 * sw_rate_1;
    }
    else
    {
      p_0 = p_0 * sw_rate_0;
      p_1 = p_1 * (1 - sw_rate_1);
    }
  }
  
  u = runif(1)[0];
  
  if(u < p_0 / (p_0 + p_1))
  {
    res = 0;
  }
  else
  {
    res = 1;
  }

  return res;
}


//Precompiled C++ functions for use in R 

// [[Rcpp::export]]
void msm_generate_observations(int N_data)
{
  vec Y(N_data);
  vec reg(N_data);
  vec state(N_data);
  double sw_rate_0 = 0.001;
  double sw_rate_1 = 0.005;
  double sigma_0 = 1;
  double sigma_1 = sqrt(10);
  double beta = sqrt(1);
  
  double u;
  
  // generate regimes of the chain
  u = runif(1)[0];
  
  if(u<1. / (1 + sw_rate_1/sw_rate_0))
  {
    reg(0) = 0;
  }
  else
  {
    reg(0) = 1;
  }

  for(int i=1;i<N_data; i++)
  {
    u = runif(1)[0];
    if(reg(i-1)==0){
      if(u<1 - sw_rate_0) reg(i) = 0;
      else reg(i) = 1;
    }
    else
    {
      if(reg(i-1)==0){
        if(u<sw_rate_1) reg(i) = 0;
        else reg(i) = 1;
      }
    }
  }
  
  // generate hidden states
  if(reg(0)==0) state(0) = sigma_0 * rnorm(1)[0];
  else state(0) = sigma_1 * rnorm(1)[0];
  
  for(int i=1;i<N_data; i++)
  {
    if(reg(i)==0)  state(i) = state(i-1)+ sigma_0 * rnorm(1)[0];
    else state(i) = state(i-1)+ sigma_1 * rnorm(1)[0];
  }
  
  // generate observations
  for (int i=0; i<(N_data); i++) {
    Y(i) = state(i)+ beta * rnorm(1)[0];
  }
  
  state.save(working_directory+"msm_hidden_states");
  reg.save(working_directory+"msm_regimes");
  Y.save(working_directory+"msm_observations");
}

// [[Rcpp::export]]
NumericVector msm_get_observations()
{
  vec Y;
  Y.load(working_directory+"msm_observations");
  return wrap(Y);
}

#endif /* markov_switching_model_hpp */
