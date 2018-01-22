#ifndef density_list_hpp
#define density_list_hpp
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <cstdlib>
#include <RcppArmadillo.h>

#include "include/Adaptive_Gibbs.hpp"
#include "examples/gaussian_target.hpp"
#include "examples/template.hpp"
#include "examples/TMVN/TMVN.hpp"
#include "examples/poisson_hierarchical_model.hpp"
#include "examples/markov_switching_model.hpp"



using namespace Rcpp;
using namespace std;
distribution_class* density_list(String name)
{
    distribution_class* null_density = nullptr;
    
    if(name == "gaussian") return new gaussian;
    else if(name == "TMVN") return new tmvn;
    else if(name == "poisson_hierarchical_model") return new poisson_hierarchical_model;
    else if(name == "markov_switching_model") return new markov_switching_model;
    else if(name == "new_density") return new new_density;
    else
    {
        err_name = "error: no such distribution :(";
        err_type = 3;
        return null_density;
    }
}


#endif /* density_list_hpp */