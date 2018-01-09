dim = 50
N = 10^5

set_working_directory("Change to the desired directory")

adaptive_chain<-AMCMC(distribution_type="TMVN", dim=dim, N=N, gibbs_sampling=1, parallel_computation=1, batch_length=100, track_adaptation=1, display_progress=0, thin=1, start_location = rep(2,dim))
