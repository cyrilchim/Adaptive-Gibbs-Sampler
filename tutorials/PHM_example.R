dim = 50
N = 10^6

set_working_directory("Change to the desired directory")

adaptive_chain<-AMCMC(distribution_type="poisson_hierarchical_model", dim=dim, N=N, gibbs_sampling=0, logdensity = 1, full_cond = 1, parallel_computation=1, batch_length=100, track_adaptation=1, display_progress=0, thin=1)

