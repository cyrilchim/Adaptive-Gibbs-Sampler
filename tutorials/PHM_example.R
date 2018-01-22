dim = 50
N = 10^6

set_working_directory("Change to the desired directory")

set.seed(1)
phm_example_1_data(N_data = 100, N_param = dim)

# Run adaptive metropolis within adaptive gibbs sampler for the model
adaptive_chain<-AMCMC(distribution_type="poisson_hierarchical_model", dimension=dim, N=N, gibbs_sampling=0, logdensity = 1, full_cond = 1, parallel_adaptation=1, batch_length=100, track_adaptation=1, display_progress=0, thin=1)

