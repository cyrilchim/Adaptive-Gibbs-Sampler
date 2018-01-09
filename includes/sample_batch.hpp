#ifndef sample_batch_hpp
#define sample_batch_hpp
// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>
#include "Adaptive_Gibbs.hpp"
#include "sample_batch/multinomial.hpp"

using namespace Rcpp;
using namespace arma;


// acceptance probability in Metropolis algorithm computed from density ratio
double acceptance_density(vec prop, vec theta,int gp, distribution_class* c_density,
                          param* alg_param)
{
    if(alg_param->c_density_flag==0)
    {
        //if(alg_param.full_cond==0){
        alg_param->proposal_log_density  = log(
            NumericVector(as<Function>(alg_param->R_density)(prop) )[0]);
        
        if(alg_param->proposal_log_density != -INFINITY)
        {
            //recompute current log-density value if necessary
            if(alg_param->reset_log_density == 1)
            {
                alg_param->current_log_density = log(
                    NumericVector( as<Function>(alg_param->R_density)(theta) )[0]);
    
                alg_param->reset_log_density = 0;
            }
            
            return  exp(alg_param->proposal_log_density - alg_param->current_log_density);
        }
        else return 0;
        //}
        /*else
         {
         b_prop  = as<Function>(alg_param.R_density)(wrap(prop), gp);
         b =  as<Function>(alg_param.R_density)(wrap(theta), gp);
         return  log(b_prop[0]) - log(b[0]);
         }*/
    }
    
    else
    {
        if(alg_param->full_cond==0)
        {
            alg_param->proposal_log_density = log(c_density->density(prop));
            
            if(alg_param->proposal_log_density != -INFINITY)
            {
                //recompute current log-density value if necessary
                if(alg_param->reset_log_density == 1)
                {
                    alg_param->current_log_density =  log(c_density->density(theta));
                    alg_param->reset_log_density = 0;
                }
                return  exp(alg_param->proposal_log_density - alg_param->current_log_density);
            }
            else return 0;
        }
        else
        {
            alg_param->proposal_log_density = log(c_density->full_cond(prop, gp));
            if(alg_param->proposal_log_density != -INFINITY)
            {
                //recompute current conditional log-density value if necessary
                alg_param->current_log_density =  log(c_density->full_cond(theta, gp));
                
                return  exp(alg_param->proposal_log_density - alg_param->current_log_density);
            }
            else return 0;
        }
        
    }
}

// acceptance probability in Metropolis algorithm computed from log-density ratio
double acceptance_logdensity(vec prop, vec theta,int gp, distribution_class* c_density,
                             param* alg_param)
{
    if(alg_param->c_density_flag==0)
    {
        //if(alg_param.full_cond==0){
        alg_param->proposal_log_density  = NumericVector(
          as<Function>(alg_param->R_density)(prop))[0];
        
        if(alg_param->proposal_log_density != -INFINITY)
        {
            //recompute current log-density value if necessary
            if(alg_param->reset_log_density == 1)
            {
                alg_param->current_log_density =  NumericVector(
                  as<Function>(alg_param->R_density)(theta))[0];
              
                alg_param->reset_log_density = 0;
            }
            
            return  exp(alg_param->proposal_log_density - alg_param->current_log_density);
        }
        else return 0;
    }
    
    else
    {
        if(alg_param->full_cond==0)
        {
            alg_param->proposal_log_density = c_density->logdensity(prop);
            
            if(alg_param->proposal_log_density != -INFINITY)
            {
                //recompute current log-density value if necessary
                if(alg_param->reset_log_density == 1)
                {
                    alg_param->current_log_density =  c_density->logdensity(theta);
                    alg_param->reset_log_density = 0;
                }
                return  exp(alg_param->proposal_log_density - alg_param->current_log_density);
            }
            else return 0;
        }
        else
        {
            alg_param->proposal_log_density = c_density->logfull_cond(prop, gp);
            if(alg_param->proposal_log_density != -INFINITY)
            {
                //recompute current conditional log-density value if necessary
                alg_param->current_log_density =  c_density->logfull_cond(theta, gp);
                
                return  exp(alg_param->proposal_log_density - alg_param->current_log_density);
            }
            else return 0;
        }
        
    }
}



//RSGS STEP (Step 1 from ARSGS paper)
// here theta is a starting location for a batch, scale - vector of scales in Metropolis-within-Gibbs algorithm
// tr - trace of the chain
// tr_aux - auxiliary variable that records samples that are used by adaptive_step in order to update the covariance S
// start - index of the variable from which the batch is sampled
// count - batch number
void sample_batch(vec* theta, vec* scale, mat *L_1, field<mat>* S_cond, unsigned int* count,
                  vector<vector<double> > *tr, vector<vector<double> > *tr_aux, int *start,
                  distribution_class* c_density, param *alg_param)
{
    int i, proposal_update;
    proposal_update = ceil( double(alg_param->batch_length) / double(alg_param->frequency_ratio) );
    
    double  gp, A, u;
    vec prop = *theta;
    vec accum(par);
    uvec updates_number(par);
    
    unsigned bl_start, bl_end;
    
    accum.zeros();
    updates_number.zeros();
    
    // gsl_ran_discrete_t *g;
    // g = gsl_ran_discrete_preproc(par, alg_param->p);
    
    if(alg_param->burn_in==1){alg_param->batch_length = burn_in;}
    
    for (int l = 1; l<=alg_param->batch_length; l++)
    {
        for (i = 0; i<par; i++)
        {
            gp = multinomial(alg_param->p); // choose a block using the sampling weights
            // gp = gsl_ran_discrete(r, g);
            // gp  = i;
            
            if((alg_param->gibbs_step[int(gp)]==0&&alg_param->ignore_gibbs_sampling==1) ||
               (alg_param->gibbs_sampling==0&&alg_param->ignore_gibbs_sampling==0))
            {
                //Metropolis adaptation if chosen
                
                //generate a proposal
                bl_start = alg_param->blocking_structure[int(gp)];
                bl_end = alg_param->blocking_structure[int(gp)+1]-1;
                
                u = runif(1)[0];
                if(u < alg_param->stabilizing_weight)
                  {
                    prop.subvec(bl_start, bl_end) = (theta->subvec(bl_start, bl_end) +
                      (*scale)(gp) * (*S_cond)(gp) * vec( rnorm(bl_end - bl_start + 1) ));
                  
                  }
                else
                  {
                  prop.subvec(bl_start, bl_end) = (theta->subvec(bl_start, bl_end)
                    + 2.38/sqrt(dim)  * vec( rnorm(bl_end - bl_start + 1) ));
                  }
                
                
                //compute logarithm of acceptance ratio
                if(alg_param->logdensity==0)
                {
                    A = acceptance_density(prop, *theta, gp, c_density, alg_param);
                }
                else
                {
                    A = acceptance_logdensity(prop, *theta, gp, c_density, alg_param);
                }
                
                if(1<A) {A = 1;}
                
                //accumulated acceptance ratio as in the AIRMCMC algorithm
                if(alg_param->adapt_proposals==1)
                {
                    accum(gp) = accum(gp) + A;
                    updates_number(gp) = updates_number(gp) + 1;
                }
                
                
                //accept-reject step
                bl_start = alg_param->blocking_structure[int(gp)];
                bl_end = alg_param->blocking_structure[int(gp)+1]-1;
                
                u =  runif(1)[0];
                
                if(u<=A)
                {
                    //accept proposal
                    theta->subvec(bl_start, bl_end)  = prop.subvec(bl_start, bl_end) ;
                    //change current value of log_density (or conditional log density)
                    if(alg_param->full_cond == 0)
                    {
                        alg_param->current_log_density = alg_param->proposal_log_density;
                    }
                    //change value of the current_location parameter
                    current_location.subvec(bl_start, bl_end) = (
                      theta->subvec(bl_start, bl_end));
                }
                else
                {
                    //reject proposal
                    prop.subvec(bl_start, bl_end)  = theta->subvec(bl_start, bl_end);
                }
                
            }
            
            
            if(alg_param->gibbs_step[int(gp)]==1)
            {
                //Gibbs adaptation if chosen
                theta->subvec(alg_param->blocking_structure[int(gp)], alg_param->blocking_structure[int(gp)+1]-1) = c_density->sample_full_cond(*theta, gp);
                alg_param->reset_log_density = 1;
                
                //update the proposal vector for the Metropolis step
                if(alg_param->gibbs_sampling==0)
                {
                    prop.subvec(alg_param->blocking_structure[int(gp)], alg_param->blocking_structure[int(gp)+1]-1)  = theta->subvec(alg_param->blocking_structure[int(gp)], alg_param->blocking_structure[int(gp)+1]-1) ;
                }
                
                //change value of the current_location parameter
                current_location.subvec(alg_param->blocking_structure[int(gp)], alg_param->blocking_structure[int(gp)+1]-1) = theta->subvec(alg_param->blocking_structure[int(gp)], alg_param->blocking_structure[int(gp)+1]-1);
                
            }
        }
        
        
        
        //adapt proposal variances (scales)
        if(alg_param->adapt_proposals==1 && alg_param->burn_in==0 && l % proposal_update == 0) {
            for(int j=0;j<par;j++){
                if(alg_param->blocking_structure(j+1) - alg_param->blocking_structure(j) == 1 && updates_number(j)>0){
                    (*scale)(j) = (*scale)(j)*exp( pow(*start+l,-0.7/(1+alg_param->beta) )  * (accum(j)/updates_number(j)-0.44));
                }
                else if(alg_param->blocking_structure(j+1) - alg_param->blocking_structure(j) > 1 &&updates_number(j)>0){
                    (*scale)(j) = (*scale)(j)*exp( pow(*start+l,-0.7/(1+alg_param->beta) )  * ( accum(j)/updates_number(j)-0.234 )  );
                }
            }
            
            accum.zeros(); updates_number.zeros();
        }
        
        
        //update the number of sampled points
        if(alg_param->burn_in==0){*start = *start + 1;}
        
        //record sampled point
        if((alg_param->adapt_proposals==1 || alg_param->adapt_weights==1 || alg_param->estimate_spectral_gap==1)&&(alg_param->burn_in==0))
        {
            for(i=0;i<dim;i++) 
            {
                (*tr_aux)[i].push_back((*theta)(i));
            }
        }
        //(*tr_aux)[i][l - 1] = (*theta)(i);
        if(alg_param->save==1 && (*start % alg_param->thin==0)) {
            for(i=0;i<dim;i++) 
            {
                (*tr)[i].push_back((*theta)(i));
            }
        }
        
        
    }
    
    //update the number of sampling weights adaptations
    if(alg_param->burn_in==0){*count = *count + 1;} 
}

#endif /* sample_batch_hpp */