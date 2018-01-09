
RcppParallel::setThreadOptions(numThreads = 2) #we will need 2 threads to run the parallel version of ARSGS

#Example 1. The same as the first example from tutorian_examples_test.R 

# set working directory to a folder, where the downloaded files are
setwd("change this to a directory where the library files are")
#setwd("change this to a desired directory")


# compile the library to a desired folder
Rcpp::sourceCpp('Adaptive_Gibbs_parallel.cpp',cacheDir = "change this to a desired existing directory")
#Rcpp::sourceCpp('Adaptive_Gibbs.cpp',cacheDir = "change this to a desired directory")


# define a directory where the output of the algorithm will be stored
directory <- "change this to a desired existing directory"
#direcotry <- "change this to a desired directory"

#set random generator seed

set.seed(1)

# check if the directory exists
if(file.exists(directory)){
  # set up the directory for the output
  set_working_directory(directory)
  # set dimenssionality of the target distribution dim, and the number of desired samples N
  dim <- 10
  N <- 10000
  # set a covariance matrix such as in Example of ARSG paper (which equals to a correlation matrix)
  set_example_covariance(dim)
  # this is how the covariance matrix look like:
  print(get_covariance())
  cat("",sep = "\n\n")
  # note that set_example_covariance(dim) and get_covariance() are defined in gaussian_target.hpp file
  # run an adaptive MCMC for the gaussian distribution with the covariance matrix as above using full dimensional Random Scan Gibbs Sampler
  adaptive_chain<-AMCMC(distribution_type = "gaussian", dim = dim, N = N, gibbs_sampling = 1, parallel_computation = 1)
}else{
  print("The directory does not exist. Please set up an existing path")
}

# let's compute the estimated correlation matrix

S <- matrix(0,nrow = N, ncol = dim)

for(i in 1:dim)
{
  v <- trace_coord(i)
  S[,i] <- v
}

S <- cor(S)
# estimated correlation matrix:
S


# estimated optimal weigts
adaptive_chain$weights


# estimated spectral gap:
adaptive_chain$sp_gap



Prob <- trace_weights()
# trace of the estimated optimal weigts:
Prob

Gap <- trace_inv_sp_gap()
# trace of the estimated 1/spectral gap:
Gap



#Example 2. Here we compare execution time of the parallel algorithm vs non-parallelised on a 500 - dimensional example. We sample 10000 points using Gibbs moves with blocking c(10, 190, 300)


if(file.exists(directory)){
  # set up the directory for the output
  set_working_directory(directory)
  # set dimenssionality of the target distribution dim, and the number of desired samples N
  dim <- 500
  N <- 10000
  set_example_covariance(dim)
  # run the adaptive MCMC. Here we start the algorithm at a point (1,..,1) with starting proposal variances of the blocks (5,5,5)
  ptm<-proc.time()
  adaptive_chain<-AMCMC(distribution_type = "gaussian", dim = dim, N = N, full_cond = 1, blocking = c( rep(1,80), rep(4,105) ), save = 0, adapt_proposals = 0, adapt_weights = 0, estimate_spectral_gap = 0, parallel_computation = 0)
  print("Time spent without adaptations: ")
  proc.time() - ptm
}else{  
  print("The directory does not exist. Please set up an existing path")
}





if(file.exists(directory)){
  # set up the directory for the output
  set_working_directory(directory)
  # set dimenssionality of the target distribution dim, and the number of desired samples N
  dim <- 500
  N <- 10000
  set_example_covariance(dim)
  # run the adaptive MCMC. Here we start the algorithm at a point (1,..,1) with starting proposal variances of the blocks (5,5,5)
  ptm<-proc.time()
  adaptive_chain<-AMCMC(distribution_type = "gaussian", dim = dim, N = N, full_cond = 1, blocking = c( rep(1,80), rep(4,105) ), save = 0, parallel_computation = 0)
  print("Time spent without parallelisation: ")
  proc.time() - ptm
}else{  
  print("The directory does not exist. Please set up an existing path")
}





if(file.exists(directory)){
  # set up the directory for the output
  set_working_directory(directory)
  # set dimenssionality of the target distribution dim, and the number of desired samples N
  dim <- 500
  N <- 10000
  set_example_covariance(dim)
  # run the adaptive MCMC. Here we start the algorithm at a point (1,..,1) with starting proposal variances of the blocks (5,5,5)
  ptm<-proc.time()
  adaptive_chain<-AMCMC(distribution_type = "gaussian", dim = dim, N = N, full_cond = 1, blocking = c( rep(1,80), rep(4,105) ), save = 0, parallel_computation = 1)
  print("Time spent with parallelisation: ")
  proc.time() - ptm
}else{  
  print("The directory does not exist. Please set up an existing path")
}



