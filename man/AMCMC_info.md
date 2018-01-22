# Description

The package provides a set of functions to run Adaptive Random Scan Gibbs (ARSG) and  Adaptive Metropolis-within-Gibbs (AMwG) algorithms for a user-defined target distribution. The adaptive algorithms are defined in the ARSG paper. The adaptations tune the sampling weights of the algorithms, and also the variances of the proposals in the AMwG algorithm.  The algorithms are implemented in C++ providing R interface. While the user can define the target density in Rm it is strongly recommended to use a C++ template provided with the package in [template](../examples/template.hpp).  Gibbs Sampling is implemented only for C++ defined densities, where the user has to explicitly define a function that samples from full conditionals. The main function AMCMC(...) produces a trace of the chain which  is saved as a set of text files in a pre-specified  folder. Please refer to appendix.txt for details for additional functions provided by the package. 
			

## Remark

It is possible to simultaneously sample and adapt the parameters as described in Section 8.3 of ARSGS paper, but only for C++ defined target distributions.


### Usage
```R
List AMCMC(string distribution_type = "gaussian", int N = 1000, nullable_func R_density=R_NilValue,
           double frac_burn_in = 10, int thin = 1, double batch_length = 100,
           int frequency_proposal_update = 100,  NumericVector start_location = NA_REAL, int dimension = 1,
           bool adapt_proposals = 1, bool adapt_weights = 1, NumericVector start_weights=NA_REAL,
           bool logdensity = 0, bool full_cond = 0,
           bool gibbs_sampling = 0, bool estimate_spectral_gap = 1, double rate_beta = 0,
           NumericVector blocking = NA_REAL,  NumericVector gibbs_step = NA_REAL,
           NumericVector start_scaling = NA_REAL, int precision = 8, bool save = 1,
           bool reweight = 0, double stabilizing_weight = 0.1, double stabilizing_coeff = NA_REAL,
           bool perturb_covariance = 0, double perturb_coeff = NA_REAL, bool track_adaptation = 0,
           bool display_progress = 1, bool parallel_adaptation = 0)
```

### Arguments

  - `distribution_type` - the target distribution with its attributes defined in C++. By default, "gaussian" target distribution is set, where the covariance matrix can be set using set_covariance(matrix) or set_example_covariance(dim) functions. See tutorial_parallel_examples.R and how_to_write_density.txt files. To sample from a user-defined distribution, set distribution_type  = "new_density"
  
  - `N` - number of samples to be produced by the adaptive chain. By default, `N = 1000`;

  - `R_density` - target density defined in R. By default it is set to a null value. If logdensity is enabled, R_density is treated as a log-density of the target distribution;

  - `frac_burn_in` - length of burn-in run. Should be a fraction of the total length N. By default, frac_burn_in =  10 percent, i.e., the default length of the burn-in is `0.1*N`;
	
  - `thin` - thinning parameter. By default, `thin = 1`;

  - `batch_length` - number of iterations between the adaptation of the sampling weights. By default, batch_length = 100. This value may be tuned by `rate_beta` parameter;

  - `frequency_proposal_update` - relative frequency of the proposal scalings (variances) adaptations. By default, frequency_ratio  = 10, meaning the proposals are updated 10 times more frequently than the weights'
	
  - `start_location` - staring location. By default, it is set to be the origin, i.e., (0,.., 0);
	
  - `dimension` - dimensionality of the target distribution. By default, `dimension = 1`;

  - `adapt_proposals`  - enable/disable proposal variances adaptation. Enabled by default;

  - `adapt_weights` - enable/disable weights adaptation. Enabled by default;
	
  - `start_weights` - starting probability weights for the Random Scan Gibbs Sampler;

  - `logdensity` - use the logarithm of the density when computing the acceptance ratio in the Metropolis part of the algorithm;
	
  - `full_cond` - whether to use manually defined full conditional densities when performing Metropolis-within-Gibbs algorithm. By default, `full_cond = 0`;

  - `gibbs_sampling` - switch between Random Scan Gibbs Sampler ( = 1) if defined in "template.hpp", and Adaptive Metropolis-within-Gibbs (= 0);

  - `estimate_spectral_gap` - enable/disable estimation of the pseudo-spectral gap defined in the ARSG paper;

  - `rate_beta` - non-negative constant beta defined in Theorem 1 of the AIRMCMC paper; after every sampling weights adaptation, updates batch_length  =  batch_length * count^ rate_beta, where count is the sampling weights adaptation number;

  - `blocking` - blocking structure for the algorithm. For example, if dim  = 10 and blocking  = c(5 ,2, 3), coordinates are grouped into three blocks: first 5 coordinates are combined into a first block, coordinates 6 and 7 into a second block, and 8, 9, and 10 into a third block. The algorithm then updates all coordinates within a block simultaneously;

  - `gibbs_step` - a vector of 0/1 of the same length as blocking. Every entry indicates wether the corresponding block is updated using its full conditional distribution ( = 1) or using a Metropolis update with a proposal  from a normal distribution. For example, say dim  = 10, blocking  = c(5 ,2, 3), and gibbs_step  = c(0, 1, 0). Then the algorithm will update the first and third blocks using a Metropolis update with normal proposals, and the second block will be updated from its full conditional distribution;
 
  - `start_scaling` - a vector of positive real numbers of the same length as blocking. Every entry defines a starting value for the standard deviation scaling of the Metropolis proposal for the corresponding block

  - `precision` - mantissa of the values written to the output files. By default, precision = 8

  - `save` - enable/disable saving the output of the algorithm

  - `reweight` - re-weight the sampling probabilities (p_1,.., p_s) as discussed in the accompanying  paper: if block i is updated using Metropolis proposal, change p_i to d_i * p_i, where d_i is a size of the block i.
  
  - `stabilizing_weight` - a real valued number between [0,1) used to stabilise the proposal: proposal = stabilizing_weight * N(current_state, stabilizing_coeff^2) + (1 - stabilising_weight) * adaptive_proposal;
  
  - `stabilizing_coeff` -  positive real number used by for generating a proposal, see `stabilizing_weight`;
  
  - `perturb_covariance` -  enable/disable perturbance of the covarince estimation: S = S + perturb_coeff * I;
  
  - `perturb_coeff` -  perturbation coefficient used by `perturb_covariance`;
  
  - `track_adaptation` - whether to print estimated values of the inverse spectrac gap and the optimal selection probabilities on the fly, after every adaptation of the selection probabilities;

  - `display_progress` - enable/disable progress bar on the screen. When enabled, might slow down the algorithm;

  - `parallel_adaptation`  - enable/disable parallelisation.


### Value

  - Returns a list with a value of estimated pseudo-spectral gap and pseudo-optimal sampling weights discussed in ARSG paper. Also, if argument `save = 1`, saves to a pre-specified directory (see set_working_directory(directory) function in appendix.txt) the trace of each coordinate, the trace of adapted proposal scaling parameters, the trace of adapted sampling weights, and the trace of the adapted values of 1 over pseudo-spectral gap.









