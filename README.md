# Adaptive-Gibbs-Sampler
Implementation of the Adaptive Gibbs and Adaptive Metropolis within Adaptive Gibbs Samplers.

## Description
The package provides a set of functions to run Adaptive Random Scan Gibbs (ARSG) and  Adaptive Metropolis-within-Gibbs (AMwG) algorithms for a user-defined target distribution. The adaptive algorithms are defined in the ARSG paper. The adaptations tune the sampling weights of the algorithms, and also the variances of the proposals in the AMwG algorithm.  The algorithms are implemented in C++ providing R interface. While the user can define the target density in Rm it is strongly recommended to use a C++ template provided with the package in "template.hpp".  Gibbs Sampling is implemented only for C++ defined densities, where the user has to explicitly define a function that samples from full conditionals. The main function AMCMC(...) produces a trace of the chain which  is saved as a set of text files in a pre-specified  folder. Please refer to appendix.txt for details for additional functions provided by the package. 

This package was tested on OS X 10.10 and Ubuntu 16.04

Table of Contents
-----------------

   * [Adaptive-Gibbs-Sampler](#adaptive-gibbs-sampler)
      * [Description](#description)
      * [Installation](#installation)
         * [Compile Adaptive_Gibbs.cpp](#compile-adaptive_gibbscpp)
      * [Using the library](#using-the-library)
         * [Gaussian example](#gaussian-example)
         * [How to write custom densities](#how-to-write-custom-densities)
	 

## Installation
Mac Users:
1. Make sure the newest version of R is installed. If RStudio is used, please update to at least version 0.10.3
2. Open terminal:
  * Run <br/>
  `xcode-select --install`<br/>
  to install xcode command line tools (this will install gcc in particular).
  *  If gfortran is not installed on your machine, run <br/>
		`curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2`<br/>
		`sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /`<br/>
    This should update your gfortran compiler (see here for details)


Ubuntu Users:
1. Make sure the last version of R is installed.
2. Make sure you have installed gcc:

	sudo apt-get install gcc

### Compile Adaptive_Gibbs.cpp

Install Rcpp, RcppArmadillo, RcppParallel, if not yet installed:
1. Open R/RStudio session and install the following packages if not yet installed:<br/>
```R
install.packages("Rcpp","RcppArmadillo", "RcppParallel")
```
2. Download/clone the repository files to a local directory, say, "~/Downloads".

3. Open R/Rstudio session use `Rcpp::sourceCpp` to compile the library:
```R
Rcpp::sourceCpp("~/Downloads/Adaptive-Gibbs-Sampler/Adaptive_Gibbs_parallel.cpp")
```

3. If you want to compile the library later to a specific folder, say, "~/Downloads/Adaptive-Gibbs-SamplerAdaptiveGibbsLib", you need to specify `cacheDir` argument:
```R
Rcpp::sourceCpp("~/Downloads/Adaptive-Gibbs-Sampler/Adaptive_Gibbs_parallel.cpp", cacheDir = "~/Downloads/Adaptive-Gibbs-Sampler/AdaptiveGibbsLib")
```
4. Use Step 3 to load precompiled library. It will be recompiled automatically if any changes are made to the source files.

## Using the library

After loading.compiling the library:
1. First, set up a folder where all the output will be stored using `set_working_directory(path)`, where `path` is the desired folder, say 
```R
set_working_directory("Users/Admin/AdaptiveGibbs/simulation_results/")
```
Notice, that you need to indicate the precise path to the directory, NOT "~/AdaptiveGibbs/simulation_results". Also, notice "/" at the end of the pass.

2. Now you can use `AMCMC(...)` and related functions. Please refer to [AMCMC_info.md](../master/man/AMCMC_info.md), [AMCMC_appendix.md](../master/AMCMC_appendix.md),  and [tutorials](../master/man/tutorials) to learn how to use `AMCMC(...)` function. 

### Gaussian example
Consider sampling from a 10-dimensional multivariate normal target distribution with a covariance matrix `S`. As an example, the user can generate a block diagonal matrix using `set_example_covariance`:
```R
dim <- 10 # dimensionality of the target 
N <- 10000 # number of desired samples
# set an example covariance matrix
set_example_covariance(dim)
# print the matrix on screen
get_covariance()
# set random seed
set.seed(42)
```
In order to sample from the distribution we need to provide at least a density function. The target distribution is described in [gaussian_target.hpp](../master/examples/gaussian_target.hpp). Sampling can be performed usinf
```C++
adaptive_chain<-AMCMC(distribution_type = "gaussian", logdensity = 1, dim = dim, N = N)
```
The function calls the coordinate-wise Adaptive Metropolis within Adaptive Gibbs Sampler for `gaussian` distribution. `logdensity` specifies that the log-density of the target distribution is used for Metropolis acceptance-rejection step. The output of the chain is written to the prespecified directory (use `get_working_directory()` to print the directory on screen). `adaptive_chain` stores the value of estimated inverse pseudo-spectral gap and the optimal sampling selection probabilities. 

Use `trace_coord(i)` to access the trace of the coordinate `i`. For example, in order to estimate the covariance matrix, one could use the following code
```R
S <- matrix(0,nrow = N, ncol = dim)
for(i in 1:dim)
  {
   v <- trace_coord(i)
   S[,i] <- v
S <- cor(S)
# estimated correlation matrix:
S
```

### How to write custom densities

Detailed guide on writing user-defined densities is described in [write_custom_density.md](../master/man/write_custom_density.md).  
One could provide an R-defined density, say,
```R
example_logdensity<-function(x)
{ 
  return( -1./2* (t(x) %*%Q_0%*%  x)[1,1] )
}
```
and call 
```R
AMCMC(R_density = example_logdensity, logdensity = 1,
                        dim = dim, N = N,  blocking = c(1,2,2,5))
```			

However, a much faster algorithm is achieved if a C++ density is supplied. A template file is provided in [template.hpp](../master/examples/template.hpp). You have to specify at least one function of the template class, say, log-density function:
```C++
double gaussian::logdensity(vec theta)
	{
	 return -1./2*dot(theta.t()*Q,theta);
	}
```
Precision matrix `Q` should be specified in the constructor. It could be read from a file. Various examples can be found [here](../master/exampes) for 
