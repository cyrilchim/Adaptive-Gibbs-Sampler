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
         * [Minimum working example](#minimum-working-example)
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

2. Now you can use `AMCMC(...)` and related functions. Please refer to [AMCMC_info.md](../master/AMCMC_info.md), [AMCMC_appendix.md](../master/AMCMC_appendix.md),  and [tutorials](../master/tutorials) to learn how to use `AMCMC(...)` function. 

### Minimum working example

### How to write custom densities

