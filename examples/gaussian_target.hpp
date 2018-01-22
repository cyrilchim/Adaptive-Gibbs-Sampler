#ifndef gaussian_target_hpp
#define gaussian_target_hpp
#define ARMA_DONT_PRINT_ERRORS
// [[Rcpp::depends(RcppArmadillo)]]
#include "includes/Adaptive_Gibbs.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;


//GAUSSIAN TARGET
class gaussian: public distribution_class{
public:
    //add some additional parameters if necessary
    mat S, I, Q, A, D; // S - covariance matrix of the target. Q - precision matrix. I - identity matrix
    field<mat> A_block;
    field<mat> sd, inv_sd;
    //add some additional functions if necessary
    gaussian();//specify constructor  if necessary
    ~gaussian(); //specify destructor if necessary
    double density(vec theta);    //target density at the point theta
    double logdensity(vec theta); //logarithm of the target density
    
    double full_cond(vec theta, int ind);  //full conditional of the block `ind`
    double logfull_cond(vec theta, int ind);  //logarithm of the full conditional of the block `ind`
    
    vec sample_full_cond(vec theta, int ind);   //sample block `ind` from its full conditional distributiont
};


gaussian::gaussian(){
    S.set_size(dim,dim); 
  
   
    
    if(S.load(working_directory+"gaussian_covariance.mat")==0){
      err_name = "error: can't initialise gaussian distribution: file not found. Use set_covariance or set_example_covariance";
      err_type = 2;
    }
    else{
      S.load(working_directory+"gaussian_covariance.mat");
      if(S.n_cols!=dim || S.n_rows!=dim)
      {
        err_name = "error: can't initialise gaussian distribution: wrong dimensionality of the input covariance matrix";
        err_type = 2;
      }
    }
    
  
    if(err_type == 0)
    {
      I.set_size(dim,dim);
      I.eye();
      D.set_size(dim, dim);
      D.zeros();
      
      Q.set_size(dim,dim);
      Q = inv(S);
      
      A.set_size(dim,dim);
      A_block.set_size(par,par);
   
      sd.set_size(par);
      inv_sd.set_size(par);
      

      
      for(int i = 0; i<par; i++){
        D.submat(blocking_structure(i),blocking_structure(i),blocking_structure(i+1)-1,blocking_structure(i+1)-1) = Q.submat(blocking_structure(i),blocking_structure(i),blocking_structure(i+1)-1,blocking_structure(i+1)-1); //update diagonal block matrix D^{-1} from the algorithm
      }
      
      A = I - inv(diagmat(D))*Q;
      for(int i=0;i<par;i++){
        sd(i) = inv(chol(D.submat(blocking_structure(i),blocking_structure(i),blocking_structure(i+1)-1,blocking_structure(i+1)-1))); 
        inv_sd(i) = chol(D.submat(blocking_structure(i),blocking_structure(i),blocking_structure(i+1)-1,blocking_structure(i+1)-1));
        }
      
        for(int i = 0; i<par; i++)
        {
          A.submat(blocking_structure(i),blocking_structure(i),blocking_structure(i+1)-1,blocking_structure(i+1)-1).zeros();
          for(int j = 0; j<par; j++)
          {
            A_block(i,j) = A.submat(blocking_structure(i),blocking_structure(j),blocking_structure(i+1)-1,blocking_structure(j+1)-1); //update diagonal block matrix D^{-1} from the algorithm
          }
        }
    }
    //cout<<"Q: "<<Q<< "inv_sd: "<<inv_sd<<"\n";
}


gaussian::~gaussian(){}



double gaussian::density(vec theta)
{
    return exp(-1./2*dot(theta.t()*Q,theta));
    
}

double gaussian::logdensity(vec theta){return -1./2*dot(theta.t()*Q,theta);}




double gaussian::full_cond(vec theta, int ind)
{
    vec mn(blocking_structure(ind+1) - blocking_structure(ind)), v(blocking_structure(ind+1)-blocking_structure(ind)); 
    mn.zeros();
   
    for (int j = 0; j<par; j++)
    {
      mn = mn + A_block(ind, j) * theta.subvec(blocking_structure(j), blocking_structure(j+1)-1);
    }

   v = inv_sd(ind) * (theta.subvec(blocking_structure(ind), blocking_structure(ind+1)-1) - mn);
   return exp(-1./2*dot(v,v));
}


double gaussian::logfull_cond(vec theta, int ind)
{

  vec mn(blocking_structure(ind+1) - blocking_structure(ind)), v(blocking_structure(ind+1)-blocking_structure(ind)); 
  mn.zeros();
  
  for (int j = 0; j<par; j++)
  {
    mn = mn + A_block(ind, j) * theta.subvec(blocking_structure(j), blocking_structure(j+1)-1);
  }
 
  v = inv_sd(ind) * (theta.subvec(blocking_structure(ind), blocking_structure(ind+1)-1) - mn);
  return -1./2*dot(v,v);
}


vec gaussian::sample_full_cond(vec theta, int ind)
{
    vec res(blocking_structure(ind+1) - blocking_structure(ind)); 
    res.zeros();
    
    for (int j = 0; j<par; j++)
    {
        res = res + A_block(ind, j) * theta.subvec(blocking_structure(j), blocking_structure(j+1)-1);
        //res = res + A(ind, j) * theta(ind);
      
    }
    res = res + sd(ind) * vec(rnorm(blocking_structure(ind+1) - blocking_structure(ind)));

    return  res;
   
}

// [[Rcpp::export]]


void set_example_covariance(int dim)
{
  mat S;
  S.set_size(dim,dim);
  S.eye();
  if(dim>=2){
    S(1,0) = S(0,1) = - 0.95;
    for(int i=2;i<=dim/2;i++) S(2*i-2,2*i-1) = S(2*i-1, 2*i-2) = - 0.95/i;
  }
  S.save(working_directory+"gaussian_covariance.mat");
}

// [[Rcpp::export]]
void set_covariance(NumericMatrix S_save)
{
    int dim = S_save.ncol();
    mat S(dim, dim);
    for(int i=0;i<dim;i++)
        for(int j=0;j<dim;j++)
            S(i,j) = S_save(i,j);
    S.save(working_directory+"gaussian_covariance.mat");
}

// [[Rcpp::export]]
NumericMatrix get_covariance()
{
    mat S;
    S.load(working_directory+"gaussian_covariance.mat");
    return wrap(S);
}

#endif /* gaussian_target_hpp */
