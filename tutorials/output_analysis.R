library(batchmeans) #estimate asymptotic variance for MCMC

set_working_directory("path to the file with the simulation results")

w <- trace_weights()

par = length(w[1,])

par = 50

asymp_var<-rep(0, par)
eff_samp_size<-rep(0,par)


Sigma = estimated_covariance()
D<-diag(Sigma)


for(i in 1:par)
{
  v <- trace_coord(i)
  eff_samp_size[i] = ess(v)
  asymp_var[i] = (bm(v)$se)^2 * length(v)/var(v)
}

# Maximum asymptotic variance amongst coordinate projection function f(x) = x_i/sqrt(Var(x_i))
max(asymp_var)


max(w[length(w[,1]),])
# Maximum potential speed up against the vanilla Ransom Scan Gibbs Sampler
1./2 * max(w[length(w[,1]),]) * dim

