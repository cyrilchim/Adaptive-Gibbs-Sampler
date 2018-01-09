// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppGSL)]]

#define ARMA_DONT_PRINT_ERRORS

#include "density_list.hpp"

#include "includes/adaptive_step.hpp"
#include "includes/sample_batch.hpp"
#include "includes/reweight.hpp"

#include <RcppParallel.h>

#include "includes/progress_bar.hpp"

// #include <RcppArmadilloExtensions/sample.h>
// #include <gsl/gsl_rng.h>
// #include <gsl/gsl_randist.h>

using namespace tbb; //support for parallel computing
using namespace Rcpp;
using namespace arma;
using namespace std;

// gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);

distribution_class* c_density;

//parallel sampling in a manner described n ARSG paper
void parallel_step(mat *S, vec* mu, int* start_old, int* start_new, vector<vector<double> >* tr_aux_tmp,  mat* Q, mat* L_1_tmp, field<mat>* S_cond_tmp, vec *scale_tmp, vec* dir, vec* z, vec* w, vec* p, double* inv_sp_gap, vector<double> *inv_sp_gap_trace, vector<vector<double> > *tr_proposals, param *alg_param_tmp, vec* theta, vec* scale, mat *L_1, field<mat>* S_cond, unsigned int* count, vector<vector<double> > *tr, vector<vector<double> > *tr_aux, int *start, param *alg_param)
{
  
  vec perturb = vec(rnorm(dim+1)); //perturbation from Step 3.1.1 of the ARSGS
  
  task_group g; //parallel computation using tbb library (included in RcppParallel)
  g.run(  [&] { adaptive_step(&perturb, S, mu, start_old, start_new, tr_aux_tmp,  Q, L_1_tmp,
                              S_cond_tmp, scale_tmp, dir, z, w, p, inv_sp_gap, inv_sp_gap_trace,
                              tr_proposals, alg_param_tmp);  } );
  g.run(  [&] { sample_batch(theta,  scale, L_1, S_cond, count, tr, tr_aux, start, c_density,
                             alg_param); } );
  g.wait();
  
}

/********************************************************************************************************************************************/

// [[Rcpp::export]]

//Adaptive MCMC agorithm exported to R
List AMCMC(string distribution_type = "gaussian", int N = 1000, nullable_func R_density=R_NilValue,
           double frac_burn_in = 10, int thin = 1, double batch_length = 100,
           int frequency_ratio = 100,  NumericVector start_location = NA_REAL, int dims = 1,
           bool adapt_proposals = 1, bool adapt_weights = 1, NumericVector start_weights=NA_REAL,
           bool c_density_flag = 1, bool logdensity = 0, bool full_cond = 0,
           bool gibbs_sampling = 0, bool estimate_spectral_gap = 1, double rate_beta = 0,
           NumericVector blocking = NA_REAL,  NumericVector gibbs_step = NA_REAL,
           NumericVector start_scales = NA_REAL, int precision = 8, bool save = 1,
           bool reweight = 1, double stabilizing_weight = 0.9,
           bool perturb_covariance = 0, bool track_adaptation = 0,
           bool display_progress = 1, bool parallel_computation = 0)
{
  
  Rcout<<"Adaptive MCMC v0.1\n"<<endl;

  //allow parallelisation for C++ densities only
  if(parallel_computation == 1 && c_density_flag == 0)
  {
    Rcout<<"Parallelization is not available for R-defined densities.\n";
    return List::create(Named("error") = 1);
  }

  if(working_directory == ""){
    Rcout<<"error: working_directory is not specified. Please use set_working_directory(...)\
      to specify the desired path\n";
    return List::create(Named("error") = 1);
  }
 

  err_type = 0;
  
  //set up burn-in and dimensionality d
  if(frac_burn_in < 0)
  {
    Rcout<<"negative value for burn-in burn-in\n";
    return List::create(Named("error") = 1);
  }
  
  
  burn_in = floor(frac_burn_in/100. * N);
  
  dim = dims;

  //number of sampled points
  int start = 0, start_old = 0, start_new = 0;

  //auxiliary variables
  int i, j, sz;
  double s;
  
 
  param alg_param;// algorithm parameters

 //check if R_density properly passed to the function
  if(c_density_flag==0&&R_density==R_NilValue)
    {
      Rcout<<"R density is not defined \n";
      return List::create(Named("error") = 1);
    }

  if(c_density_flag==0&&full_cond==1)
    {
      Rcout<<"this vesion does not support full conditionals for R_density  \n";
      return List::create(Named("error") = 1);
    }

  //SET UP blocking_structure VARIABLE
  
  NumericVector blocking_structure_tmp = NULL;

  blocking_structure_tmp.push_back(0);
  
  i=0;
  for(int j=0;j<blocking.size();j++)
  {
    i = i+blocking(j);
    blocking_structure_tmp.push_back(i);
  }
  

 if(!(NumericVector::is_na(blocking[0]))&&*prev(blocking_structure_tmp.end())!=dim)
   {
      Rcout<<"error: number of blocks does not match blocking structure\n";
      return List::create(Named("error") = 1);
    }

  if(NumericVector::is_na(blocking[0]))
  {
    //set up number of blocks
    par = dim;
    
    alg_param.blocking_structure.resize(par+1);
    for(i=0; i<=par; i++) alg_param.blocking_structure(i) = i;
  }
  else 
  {
    //set up number of blocks
    par = blocking_structure_tmp.size() - 1;
    
    alg_param.blocking_structure.resize(par+1);
    for(i=0;i<=par;i++)
    {
      alg_param.blocking_structure(i) = (unsigned int)( blocking_structure_tmp(i) );
    }
  }
  
  
  if(prod(sort(alg_param.blocking_structure)==alg_param.blocking_structure)==0 ||
     alg_param.blocking_structure(0)!=0)
    {
      Rcout<<"error: wrong blocking structure\n"; return List::create(Named("error") = 1);
    }
  
  alg_param.gibbs_step.resize(par);
  blocking_structure.set_size(par+1);   
  blocking_structure = alg_param.blocking_structure; //define a blocking structure of the algorithm


  
  //SET UP alg_param.gibbs_step
  
  if(!(NumericVector::is_na(gibbs_step[0]))&&gibbs_step.size()!=par)
    {
      Rcout<<"error: number of blocks does not match blocking structure\n";
      return List::create(Named("error") = 2);
    }
  
  if(NumericVector::is_na(gibbs_step[0]))
  {
    if (gibbs_sampling==0){
      alg_param.gibbs_step.zeros();
    }
    else
    {
      alg_param.gibbs_step.ones();
    }
    alg_param.ignore_gibbs_sampling = 0; //do not ignore gibbs_sampling parameter
  }
  else
  {
    alg_param.ignore_gibbs_sampling = 1;// ignore gibbs_sampling parameter
    for(i=0;i<par;i++)
      alg_param.gibbs_step(i) = (unsigned int)(gibbs_step(i));
  }

  //SET UP lowerb, m, am
  lowerb = 1. / (par * par);//lower bound
  m_iter = sqrt(par) * 50;
  am = log(m_iter) / (m_iter);

  //SET UP SAMPLING WEIGHTS
  vec p(par);

  if(!(NumericVector::is_na(start_weights[0]))&&start_weights.size()!=par)
    {
      Rcout<<"error: dimension of start_weights does not match dim\n";
      return List::create(Named("error") = 1);
    }
  
  if(NumericVector::is_na(start_weights[0])) //check if the passed probaility vector is null
  {
    for(i=0;i<par;i++) {p(i) = 1./par;}
  }
  else
      {
      s=0; for(i=0;i<par;i++) 
        {
        if(start_weights[i]<lowerb)
          {
            Rcout<<"error: some of the coordinates of start_weights are too small or negative\n";
            return List::create(Named("error") = 1);
          }
        s = s+start_weights[i];
        }
      for(i=0;i<par;i++) {p(i) = start_weights[i]/s;}
      }
  
  vec w(par+1);//auxiliary weights
    for(i=0;i<par;i++) w(i) = par/(1.+par)*p(i);
    w(par) = 1./(par+1);
    
    
   
    //SET UP THE REST OF THE VARIABLES  
    vec mu(dim); // mean value vector
    mu.zeros();
    vec theta(dim); //parameters for which MCMC is run
    vec scale(par); //proposal variance
    
    if(rate_beta<0)
      {
        Rcout<<"error: rate_beta has to be non-negative\n";
        return List::create(Named("error") = 1);
      }
    
    if(stabilizing_weight<0)
      {
        Rcout<<"error: `stabilizing_weight` has to be non-negative\n";
        return List::create(Named("error") = 1);
      }
    
    if(!(NumericVector::is_na(start_scales[0]))&&start_scales.size()!=par)
      {
        Rcout<<"error: dimension of start_scales does not match par\n";
        return List::create(Named("error") = 1);
      }
    
    if(NumericVector::is_na(start_scales[0])){
      scale.ones();
    }
    else{
      scale = start_scales;
      if(prod(scale>0)==0)
        {
          Rcout<<"error: proposal variances are improperly defined\n";
          return List::create(Named("error") = 1);
        }
    }
   
    mat S(dim, dim); //estimated covariance matrix
    
    mat Q(dim, dim); // estimated inverse covariance matrix
    
    mat L_1(dim, dim); //matrix L from the ASRGS algorithm computed at (1,..,1)
    L_1.eye();
    
    field<mat> S_cond; //covariance matrices of the proposals in Metropolis step of the ARSGS
    S_cond.set_size(par);
    
    vec grad_dir(par); //supergradient direction d^i from the ARSGS 
    grad_dir.zeros();
    
    vec z = vec(rnorm(dim+1)); //guess of the directional derivative

    vector<vector<double> > tr, tr_aux, tr_aux_tmp; 
    // tr - trace of the chain; tr_aux - auxiliary vectro, part of trace used to re-estimate
    // the covariance matrix; tr_aux_tmp - copy of tr_aux
    
    vector<double> trace_of_direction; //trace of a direction
    for (i = 0 ; i < dim ; i++){
      tr.push_back(trace_of_direction);
      tr_aux.push_back(trace_of_direction);
      tr_aux_tmp.push_back(trace_of_direction);
    }
     
    
    vector<vector<double> > tr_weights; //trace of the estimated optimal weights
    for (i = 0 ; i < par ; i++)
      tr_weights.push_back(trace_of_direction);
    
    vector<vector<double> > tr_proposals; //trace of the estimated optimal proposal variances
    for (i = 0 ; i < par ; i++)
      tr_proposals.push_back(trace_of_direction);
    
    vector<double> inv_sp_gap_trace; //trace of the estimated pseudo-spectral gap
    
    double inv_sp_gap = -1; //estimated pseudo-spectral gap

    //SET UP STARTING LOCATION
    if(NumericVector::is_na(start_location[0]))
      {
        theta.zeros(); //if starting location is not specified, start from the origin
      }
    else if(!(NumericVector::is_na(start_location[0])) && (start_location.size()!=dim))
      {
        Rcout<<"error: dimension of start_location does not match dim \n";
        return List::create(Named("error") = 1);
      }
    else{ 
      for(i=0;i<dim;i++) theta(i) = start_location[i];
    }
    
    current_location.set_size(dim);
    current_location = theta;

    //************************************ 
    
    //SET DISTRIBUTION TYPE
    
    alg_param.c_density_flag = c_density_flag;
    
    if(alg_param.c_density_flag == 1)
    {
      c_density = density_list(distribution_type);
      
      if(err_type!=0)
      {
        return List::create(Named(err_name) = err_type);
      }
    }

    alg_param.adapt_weights = adapt_weights;
    alg_param.adapt_proposals = adapt_proposals;
    if(c_density_flag==0) alg_param.R_density = R_density;
    else alg_param.R_density = R_NilValue;
    alg_param.batch_length = burn_in;
    alg_param.full_cond = full_cond;
    alg_param.logdensity = logdensity;
    alg_param.gibbs_sampling = gibbs_sampling;
    alg_param.estimate_spectral_gap = estimate_spectral_gap;
    alg_param.thin = thin;
    alg_param.beta = rate_beta;
    alg_param.frequency_ratio = frequency_ratio;
    alg_param.perturb_covariance = perturb_covariance;
    alg_param.stabilizing_weight = stabilizing_weight;
  
    //set up current density value

    if(alg_param.c_density_flag==0)
    {
      //if(alg_param.full_cond==0){
      alg_param.current_log_density  = log(NumericVector(
        as<Function>(alg_param.R_density)(theta))[0]);
    }
    else
    {
      if(alg_param.full_cond == 0)
      {
        if(alg_param.logdensity == 1){
          alg_param.current_log_density  = c_density->logdensity(theta);
        }
        else{alg_param.current_log_density  = log(c_density->density(theta));}
      }
        
    }
    
    //set up S_cond
    for(i=0; i<par; i++){
      (S_cond)(i) = inv(L_1.submat(alg_param.blocking_structure[i], alg_param.blocking_structure[i],
                                   alg_param.blocking_structure[i+1]-1,
                                   alg_param.blocking_structure[i+1]-1)).t();
    }

    vec perturb = vec(dim+1);//perturb vector for adaptive_step function
    
    
    alg_param.p = new double[par];
    
    // gsl_rng_set(r, 1);
    
    for(i = 0; i<par; i++) alg_param.p[i] = p(i);
    
    
    //set the estimated covariance matrix to zero
    S.zeros();
    
//************************************************************************************

    //BURN-IN
    Rcout<<"Burning..."<<endl;
    
    //number of the sampling weights adaptations
    unsigned int count = 1;

    sample_batch(&theta,  &scale, &L_1, &S_cond,  &count, &tr, &tr_aux, &start, c_density,
                 &alg_param); 
    checkUserInterrupt(); //check if interruption was called
    
  //ADAPTIVE ALGORITHM
  
  alg_param.save = save; //write to the file
  alg_param.burn_in = 0; //disable burn_in
  alg_param.batch_length = max(min(int(batch_length), N), 1);\
  
  mat L_1_tmp(dim, dim);//copy of the matrix L from the algorithm
  field<mat> S_cond_tmp; //copy of S_cond
  S_cond_tmp = S_cond;
  
  
  vec scale_tmp(par);// copy of w, p, scale
  param alg_param_tmp;//copy of alg_param
  
  Rcout<<"Sampling..."<<endl;
  
  // Initial mean estimate should be equal to the current location
  mu = current_location;
   
  // Sample first batch of coordinates
  sample_batch(&theta,  &scale, &L_1, &S_cond, &count, &tr, &tr_aux, &start, c_density, &alg_param);
  checkUserInterrupt(); //check if interruption was called to terminate R session
  
  start_new = start;
  

  
  for(i = 0; i<dim; i++) 
  {
    sz = tr_aux[i].size();

    for (j = 0; j < sz; j++) 
    {
      tr_aux_tmp[i].push_back(tr_aux[i].back());//get last element
      tr_aux[i].pop_back(); //erase last element

    }
}

  L_1_tmp = L_1;
  scale_tmp = scale;
  alg_param_tmp = alg_param;
  
  //AIRMCMC step
  alg_param.batch_length = max(int(batch_length * floor(pow(double(count), alg_param.beta))), 1);
  //check that the number of samples will not exceed N
  alg_param.batch_length = min(alg_param.batch_length, N - start);
  
  //RUN THE ADAPTIVE MCMC
  
  while(start+alg_param.batch_length<=N && alg_param.batch_length>0)
{
    
    checkUserInterrupt(); //check if interruption was called to terminate R session
    
    if(track_adaptation==1)
    {
      Rcout<<count - 2<<" inv_sp_gap: "<<inv_sp_gap<<" max(p): "<<max(p)<<" min(p): "<<min(p)<<"\n";
    }
    
    if(parallel_computation == 1)
      {
        //step of the ARSGS. Adaptation is performed in parallel to sampling
        parallel_step(&S, &mu, &start_old, &start_new, &tr_aux_tmp,  &Q,
                      &L_1_tmp, &S_cond_tmp, &scale_tmp, &grad_dir, &z, &w, &p, &inv_sp_gap,
                      &inv_sp_gap_trace, &tr_proposals, &alg_param_tmp, &theta,  &scale, &L_1,
                      &S_cond, &count, &tr, &tr_aux, &start, &alg_param);
      }
      else
        {
        //step of the ARSGS
          perturb = vec(rnorm(dim+1)); //perturbation from Step 3.1.1 of the ARSGS
          adaptive_step(&perturb, &S, &mu, &start_old, &start_new, &tr_aux_tmp, 
                        &Q, &L_1_tmp, &S_cond_tmp, &scale_tmp, &grad_dir, &z, &w, &p,
                        &inv_sp_gap, &inv_sp_gap_trace, &tr_proposals, &alg_param_tmp);

          sample_batch(&theta,  &scale, &L_1, &S_cond, &count, &tr, &tr_aux, &start,
                       c_density, &alg_param);
        }

    //display progress bar
    if(display_progress == 1)
    {
      print_bar(bar_width, start, N); 
    }
      

    if(alg_param.adapt_weights == 1)
      {
      
        for(i = 0; i<par; i++)
          {
            alg_param.p[i] = p[i];
          }
        //reweight probability weights for Metropolis update
        if(reweight == 1) {reweight_p(&alg_param); }
        
      }
    

      if(alg_param.save == 1)
      {
        //record the trace of the sampling weights
        for(i=0;i<par;i++)
        {
          (tr_weights)[i].push_back(alg_param.p[i]);
        }
      }
      
    
    //set tr_aux_tmp to be the samples generated during the last call of sample_batch
    for(i = 0; i<dim; i++) 
    {
      sz = tr_aux[i].size();

      for (j = 0; j < sz; j++) 
      {
        tr_aux_tmp[i].push_back(tr_aux[i].back());//get last element
        tr_aux[i].pop_back(); //erase last element

      }
    }
    
    
   start_new = start;
   L_1 = L_1_tmp;
   S_cond = S_cond_tmp;

   scale_tmp = scale;
   
   
   alg_param_tmp = alg_param;
   

   
   
   //AIRMCMC step
   alg_param.batch_length = max(int(batch_length * floor(pow(double(count), alg_param.beta))), 1);
   
   //check that the number of samples will not exceed N
   alg_param.batch_length = min(alg_param.batch_length, N - start);
   
   
}

// END OF THE ADAPTIVE MCMC
  
//SAVE THE OUTPUT TO working_directory 

 if(alg_param.save==1)
{
  
  
  vector<double>::iterator it; // Iterator for tr[i]
  //write a trace of the chain to a file
  ofstream trace;
  ostringstream convert;
  string file_n;
  Rcout<<"writing the chain output to "<<working_directory<<"\n";
  for(i=0;i<dim;i++)
    {
      convert << i+1;
      file_n = convert.str();
      convert.str("");
      convert.clear();
      trace.open(working_directory+"trace_of_coordinate_"+file_n+".txt");
      if(!trace.is_open())
        {
          return List::create(Named("error: could not a open file for\
                                      saving the trace of the chain") = 1);
        }
      for (it = (tr[i]).begin(); it != (tr[i]).end() ; it++) {
        trace << std::fixed << std::setprecision(precision)<< *it << endl;
        }
      trace.close();
    }
    
  //write the estimated covariance to a file 
  S.save(working_directory+"estimated_covariance.csv", csv_ascii);  

  //write weights trace to a file
  trace.open(working_directory+"weights"+".txt");
  if(!trace.is_open())
    {
      return List::create(Named("error: could not a open file for\
                                  saving the trace of the chain") = 1);
    }
  for(i=0;i<tr_weights[0].size();i++)
    {
      for (int j = 0; j < par; j++) {
        trace << std::fixed << std::setprecision(precision)<< tr_weights[j][i] << " ";
      }
      trace<<endl;
    }
  trace.close();
  //write spectral gap trace to a file
  trace.open(working_directory+"inv_spectral_gap"+".txt");
  if(!trace.is_open())
    {
      Rcout<<"could not a open file for saving the trace of the spectral gap \n;";
      return List::create(Named("error") = 1);
    }
  for(i=0;i<inv_sp_gap_trace.size();i++)
    {
      trace << std::fixed << std::setprecision(precision)<< inv_sp_gap_trace[i] << endl;
    }
  trace.close();
  
  //write proposal variances trace to a file
  trace.open(working_directory+"proposals"+".txt");
  if(!trace.is_open())
    {
      Rcout<<"could not a open file for saving the trace of the proposals \n;";
      return List::create(Named("error") = 1);
    }
  for(i=0;i<tr_proposals[0].size();i++)
    {
      for (int j = 0; j < par; j++) {
        trace << std::fixed << std::setprecision(precision)<< tr_proposals[j][i] << " ";
      }
      trace<<endl;
    }
  trace.close();
    
}

//EMPTY THE ALLOCATED MEMORY
if(alg_param.c_density_flag == 1){
  delete c_density; //delete the target distribution object (calls the destructor of the class)
}

vector< vector<double> >().swap(tr_weights);
vector< vector<double> >().swap(tr_proposals);
vector< vector<double> >().swap(tr);
vector< vector<double> >().swap(tr_aux);
vector< vector<double> >().swap(tr_aux_tmp);
vector< double >().swap(inv_sp_gap_trace);
  
NumericVector r_p(par); //return value of the weights p
for(i=0;i<par; i++){r_p[i] = alg_param.p[i];}  
  

//RETURN THE ESTIMATED OPTIMAL SAMPLING WEIGHTS AND THE PSEUDO-SPECTRAL GAP


if(alg_param.adapt_weights==1||alg_param.estimate_spectral_gap==1)
  {
    return List::create(Named("sp_gap") = 1/inv_sp_gap, Named("weights") =r_p);
  }
  else{
    return List::create(Named("weights") =r_p);
  }
 
}

// You include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
set.seed(1)



trace_coord<-function(i)
{
  work_dir = get_working_directory()
  v = read.table(paste(work_dir,"trace_of_coordinate_",i, ".txt",sep=""))
  return(v[,1])
}


trace_weights<-function()
{
  work_dir = get_working_directory()
  v = read.table(paste(work_dir,"weights",".txt",sep=""))
  return(as.matrix(v))
}

trace_proposals<-function()
{
  work_dir = get_working_directory()
  v = read.table(paste(work_dir,"proposals",".txt",sep=""))
  return(as.matrix(v))
}

trace_inv_sp_gap<-function()
{
  work_dir = get_working_directory()
  v = read.table(paste(work_dir,"inv_spectral_gap",".txt",sep=""))
  return(v[,1])
}

estimated_covariance<-function()
{
  work_dir = get_working_directory()
  S = read.table(paste(work_dir,"estimated_covariance",".csv",sep=""),sep = ",")
  return(as.matrix(S))
}


# Pointer implementation for R.
# Allows passing and changing additional parameters for R-defined target densities 
 Rpointer=function(input)
{  
  Rparameters=new.env(parent=globalenv())  
  Rparameters$value=input  
  class(Rparameters)='Rpointer'
  return(Rparameters)  
}  
*/
