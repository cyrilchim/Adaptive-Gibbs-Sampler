
# set working directory to a folder, where the downloaded files are
setwd("change this to a directory where the library files are")


# compile the library
Rcpp::sourceCpp("Adaptive_Gibbs.cpp")


# define a directory where the output of the algorithm will be stored
directory <- "change this to a desired existing directory"

# set random generator seed
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
  adaptive_chain<-AMCMC(distribution_type = "gaussian", dim = dim, N = N, gibbs_sampling = 1, display_progress = 0)
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
Prob[length(Prob[,1]),]

Gap <- trace_inv_sp_gap()
# trace of the estimated 1/spectral gap:
Gap



# Solve the same problem using mixture of Adaptive Metropolis-within-Gibbs and Adaptive RSGS with a blocking structure (2,4,4): 
# we update first block using Gibbs Sampler and the other blocks using Metropolis-within-Gibbs 
if(file.exists(directory)){
  # set up the directory for the output
  set_working_directory(directory)
  # set dimenssionality of the target distribution dim, and the number of desired samples N
  dim <- 10
  N <- 10000
  # the covariance matrix is the same:
  print(get_covariance())
  cat("",sep = "\n\n")
  # run the adaptive MCMC. Here we start the algorithm at a point (1,..,1) with starting proposal variances of the blocks (5,5,5)
  adaptive_chain<-AMCMC(distribution_type = "gaussian", dim = dim, N = N, full_cond = 1,
                        blocking = c(2,4,4), gibbs_step = c(1,0,0), rate_beta = 0.5,
                        start_scales = c(5,5,5), start_location = rep(1,10))
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
Prob[length(Prob[,1]),]

Gap <- trace_inv_sp_gap()
# trace of the estimated 1/spectral gap:
Gap


#Note, since the Gibbs Sampling move is performed for the first block, the first coordinate is ignored and is equal to 1 
Scaling <- trace_proposals()
# trace of the estimated optimal proposal standard deviation:
Scaling


#R-defined density example. Let the target density be as above
# R-defined target density should be a function of location only
# Use Rpointer if additional parameters need to be changed during the algorithm run 


# precision matrix of the target distribution
Q_0 = solve(get_covariance())


# a number of calls of the target density  
call_number = Rpointer(0)

example_logdensity<-function(x)
{ 
  call_number$value = call_number$value + 1 #change the value of the number of calls
  return( -1./2* (t(x) %*%Q_0%*%  x)[1,1] )
}


if(file.exists(directory)){
  # set up the directory for the output
  set_working_directory(directory)
  # set dimenssionality of the target distribution dim, and the number of desired samples N
  dim <- 10
  N <- 10000
  # run the Adaptive MwG algorithm for the R-defined log-density
  adaptive_chain<-AMCMC(R_density = example_logdensity, logdensity = 1,
                        dim = dim, N = N,  blocking = c(1,2,2,5))
}else{
  print("The directory does not exist. Please set up an existing path")
}

# the number of calls of the R-density function 
call_number$value

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

Prob <- trace_weights()
# trace of the estimated optimal weigts:
Prob[length(Prob[,1]),]

Gap <- trace_inv_sp_gap()
# trace of the estimated 1/spectral gap:
Gap


#Note, since the Gibbs Sampling move is performed for the first block, the first coordinate is ignored and is equal to 1 
Scaling <- trace_proposals()
# trace of the estimated optimal proposal standard deviation:
Scaling




# estimated optimal weigts
adaptive_chain$weights

# estimated spectral gap:
adaptive_chain$sp_gap



