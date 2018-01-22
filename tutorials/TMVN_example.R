dim = 50
N = 10^6

set_working_directory("Change to the desired directory")

set.seed(1)
tmvn_set_example_bounds(dim)
tmvn_generate_covariance(dim)

# Run adaptive gibbs sampler for the model
adaptive_chain<-AMCMC(distribution_type="TMVN", dimension=dim, N=N, gibbs_sampling=1, parallel_adaptation = 1, batch_length=100, track_adaptation=1, display_progress=0, thin=1, start_location = rep(2,dim))
