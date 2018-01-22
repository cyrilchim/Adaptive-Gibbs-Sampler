## Description
This is an appendix to AMCMCinfo.txt. Here we describe additional functions that are provided with [Adaptive_Gibbs.cpp](../Adaptive_Gibbs.cpp): 

```R
set_working_directory(String directory) 

String get_working_directory()

print_working_directory()

NumericVector trace_coord(int coord)

NumericMatrix trace_weights()

NumericMatrix trace_proposals()

NumericVector  trace_inv_sp_gap()

NumericMatrix estimated_covariance()

env Rpointer(input)
```


## Arguments 

  - `directory` - string, specifying a path where the output of the chain should be stored. It is recommended to specify the directory, e.g., directory = "~/AdaptiveGibbs/simulation_results". By default, directory = "", meaning that all the output is written to a root folder;
  - `coord`  - integer specifying a coordinate for which its' trace should be extracted;
  - `input` - virtually any R object;
	

## Details:
  - `set_working_directory(directory)` - assigns directory  to the paths where all the output is stored, which we call working_directory;
  - `get_working_directory()` - returns a string of the working_directory;
  - `print_working_directory()` - prints the working_directory;
  - `trace_coord(coord)` - returns a vector of a trace of a coordinate `coord`;
  - `trace_weights()` - returns a matrix each row of which is a trace of adapted probability weights as in the Adaptive Random Scan Gibbs Sampler from the ARSG paper;
  - `trace_proposals()` -  returns a matrix each row of which is a trace of adapted proposal coefficients. Here the coefficients are tuned in a way that retains the average acceptance ratio around 0.44 for one-dimensional updates, and around 0.234 for multi-dimensional updates;
  - `trace_inv_sp_gap()` - returns a vector of a trace of the estimated value for 1/pseudo-spectral gap discussed in the ARSG paper;
  - `estimated_covariance()` - returns estimated covariance matrix;
  - `Rpointer(input)` - creates an analogue of C-pointer to input. The returned object has $value attribute which returns the value of Rpointer. This function should be used  if R-defined density has external parameters that need to be changed during the algorithm run. Refer to the tutorial_examples.R for an example of usage.
