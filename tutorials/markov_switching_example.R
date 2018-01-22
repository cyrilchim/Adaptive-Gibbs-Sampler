N = 10^8
dim = 200

set_working_directory("Change to the desired directory")

set.seed(1)

# Run adaptive gibbs sampler for the model
adaptive_chain<-AMCMC(distribution_type="markov_switching_model", dimension=dim, N=N, gibbs_sampling=1, parallel_adaptation=1, batch_length=1000, frac_burn_in=1, perturb_covariance=1, track_adaptation=1, display_progress=0, thin=40)

