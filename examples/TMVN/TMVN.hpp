#ifndef TMVN_hpp
#define TMVN_hpp
#define ARMA_DONT_PRINT_ERRORS
// [[Rcpp::depends(RcppArmadillo)]]

#include "includes/Adaptive_Gibbs.hpp"
#include "rtnorm.hpp"
#include "truncated_normal.hpp"


using namespace Rcpp;
using namespace arma;
using namespace std;



//FEEL FREE TO USE BUT DO NOT CHANGE THE GLOBAL VARIABLES IN THIS FILE (THE GLOBAL VARIABLES ARE LISTED IN Adaptive_Gibbs.hpp)

class tmvn: public distribution_class{
private:
  mat S_0; // covariance matrix of the normal distribution before truncation
  mat Q_0; // precision matrix of the normal distribution before truncation
  vec D_0; // diagonal of `Q_0`
  mat I; // identity matrix of the same dimensions as `S_0`
  mat A; // `A = I - inv(diagmat(D0))*Q`
  int seed_tr_norm = 1;
  
  vec lower_trunc, upper_trunc; // lower and upper truncation bounds
  
public:
  //add some additional functions/parameters if necessary
  tmvn();//specify constructor  if necessary
  ~tmvn(); //specify destructor if necessary
  double density(vec theta);    //target density at the point theta 
  double logdensity(vec theta); //logarithm of the target density
    
  double full_cond(vec theta, int gp);  //full conditional of the block gp
  double logfull_cond(vec theta, int gp);  //logarithm of the full conditional of the block gp
  vec sample_full_cond(vec theta, int gp);  //sample block gp from its full conditional distributiont
  double truncated_normal_ab_sample(double mu, double sigma, double a, double b);
    
};


tmvn::tmvn()
{
  S_0.set_size(dim, dim);
  if(S_0.load(working_directory+"example_covariance.mat")==false || lower_trunc.load(working_directory+"lower_truncation_bound.vec")==false || upper_trunc.load(working_directory+"upper_truncation_bound.vec")==false)
    {
      err_name = "error: can't initialise gaussian distribution: file not found";
      err_type = 2;
    }
  else{
    S_0.load(working_directory+"example_covariance.mat");
    lower_trunc.load(working_directory+"lower_truncation_bound.vec");
    upper_trunc.load(working_directory+"upper_truncation_bound.vec");
    
    if(S_0.n_cols!=dim || S_0.n_rows!=dim || lower_trunc.n_elem!=dim || upper_trunc.n_elem!=dim )
    {
      err_name = "error: can't initialise gaussian distribution: wrong dimensionality or the input";
      err_type = 2;
    }
  }
  
  //Rcout<<"lower: "<<lower_trunc<<"\n";
  //Rcout<<"upper: "<<upper_trunc<<"\n";
  
  Q_0.set_size(dim, dim);
  Q_0 = inv(S_0);
  
  D_0.set_size(dim);
  D_0 = Q_0.diag();
  
  I.set_size(dim, dim);
  I.eye();
  
  A.set_size(dim, dim);
  A = I - inv(diagmat(D_0))*Q_0;
  
}

tmvn::~tmvn()
{ 
  
} 



double tmvn::density(vec theta)
{
    return 0;
}

double tmvn::logdensity(vec theta)
{
  return 0;
}

double tmvn::full_cond(vec theta, int gp)
{
  return 0;
}


double tmvn::logfull_cond(vec theta, int gp)
{
  return 0;
}


vec tmvn::sample_full_cond(vec theta, int gp)
{
  vec res(1);
  res.zeros();
  
  double var;
  
  double s = 0;
  
  //Rcout<<"A: "<<A<<"\n\n";
    
  for (int j = 0; j<dim; j++)
  {
    s = s + A(int(gp), j) * theta(j);
  }
  
  // std::pair<double, double> trunc_sample;
    double trunc_sample;
  
  var = 1.0 / sqrt(D_0(gp)); //conditional variance

  // sampling from truncated normal distribution
  //trunc_sample = rtnorm((lower_trunc(int(gp)) - s) / var, (upper_trunc(int(gp)) - s) / var);
    seed_tr_norm++;
    trunc_sample = truncated_normal_ab_sample(0., 1., (lower_trunc(int(gp)) - s) / var,
                                              (upper_trunc(int(gp)) - s) / var);
  
  //res(0) = s + var * trunc_sample.first; // extract the sampled variable
    res(0) = s + var * trunc_sample;
  return res;
}

double tmvn::truncated_normal_ab_sample(double mu, double sigma, double a, double b)

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_AB_SAMPLE samples the truncated Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double A, B, the lower and upper truncation limits.
//
//    Input/output, int &SEED, a seed for the random number
//    generator.
//
//    Output, double TRUNCATED_NORMAL_AB_SAMPLE, a sample of the PDF.
//
{
    double alpha;
    double alpha_cdf;
    double beta;
    double beta_cdf;
    double u;
    double x;
    double xi;
    double xi_cdf;
    
    alpha = ( a - mu ) / sigma;
    beta = ( b - mu ) / sigma;
    
    alpha_cdf = normal_01_cdf ( alpha );
    beta_cdf = normal_01_cdf ( beta );
    
    u = runif(1)[0];
    xi_cdf = alpha_cdf + u * ( beta_cdf - alpha_cdf );
    xi = normal_01_cdf_inv ( xi_cdf );
    
    x = mu + sigma * xi;
    
    return x;
}


// Precompiled C++ functions for use in R


// [[Rcpp::export]]
void tmvn_generate_covariance(int dim)
{
  mat S(dim, dim);
  S.eye();
  vec v(dim);
  v.zeros();
  int i;
  
  S = 0.01*S;
  
  for (i=0; i<dim; i++)
  {
    v(i) =  rbeta(1, 0.1, 0.2)[0];
    // v(i) =  1./log(2 + i) * rbeta(1, 0.1, 0.2)[0];
  }
    
  S =S + v * v.t();
  
  vec D1(dim);
  D1  = S.diag();
  D1 = sqrt(1./D1);
  
  S = diagmat(D1)*S*diagmat(D1);
  S.save(working_directory+"example_covariance.mat");
}


// [[Rcpp::export]]
void tmvn_set_example_bounds(int dim)
{
  vec upper(dim), lower(dim);
  upper.ones();
  lower.ones();
  
  upper = 3 * upper;

  lower.save(working_directory+"lower_truncation_bound.vec");
  upper.save(working_directory+"upper_truncation_bound.vec");
}

// [[Rcpp::export]]
NumericVector tmvn_get_generated_covariance()
{
  mat S;
  S.load(working_directory+"example_covariance.mat");
  return wrap(S);
}


#endif /* TMVN_hpp */
