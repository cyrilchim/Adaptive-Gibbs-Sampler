#ifndef template_hpp
#define template_hpp
#define ARMA_DONT_PRINT_ERRORS
// [[Rcpp::depends(RcppArmadillo)]]
#include "include/Adaptive_Gibbs.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;



//FEEL FREE TO USE BUT DO NOT CHANGE THE GLOBAL VARIABLES IN THIS FILE (THE GLOBAL VARIABLES ARE LISTED IN Adaptive_Gibbs.hpp)

class new_density: public distribution_class{
public:
    //add some additional functions/parameters if necessary
    new_density();//specify constructor  if necessary
    ~new_density(); //specify destructor if necessary
    double density(vec theta);    //target density at the point theta 
    double logdensity(vec theta); //logarithm of the target density
    
    double full_cond(vec theta, int gp);  //full conditional of the block gp
    double logfull_cond(vec theta, int gp);  //logarithm of the full conditional of the block gp
    vec sample_full_cond(vec theta, int gp);  //sample block gp from its full conditional distributiont
};


new_density::new_density()
{
  
}

new_density::~new_density()
{ 
  
} 



double new_density::density(vec theta)
{
  return 0;
}

double new_density::logdensity(vec theta)
{
  return 0;
}

double new_density::full_cond(vec theta, int gp)
{
  return 0;
}


double new_density::logfull_cond(vec theta, int gp)
{
  return 0;
}


//should return a vector of the same dimension as the block gp
vec new_density::sample_full_cond(vec theta, int gp)
{
  vec res;
  return res;
}


//Precompiled C++ functions for use in R 

// [[Rcpp::export]]
double some_function()
{
  return 0;
}

#endif /* template_hpp */
