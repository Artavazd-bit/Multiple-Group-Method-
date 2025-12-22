library(lavaan)
library(doParallel)
library(foreach)
library(dplyr)
library(matlib)

source("setupsmall.R")

nkernel <- 8
nobs <- c(25, 50, 100, 200, 1000)
simrunstotal <- 1000
# normal: no skewness, no kurtosis, nonnormal: specified skewness and kurtosis 
# values 
disttable <- data.frame(name = c("normal"))
disttable$skewness <- list(NULL)
disttable$kurtosis <- list(NULL)
######################### Monte Carlo Simulation ###############################
# initialising cluster 
cl <- parallel::makeCluster(nkernel, outfile = "errorcluster.txt")
doParallel::registerDoParallel(cl)

simresults <- foreach(jj = 1:nrow(simModels), .packages = c("lavaan", "foreach", "dplyr"), .combine = "rbind") %:%
  foreach(n = nobs, .combine = "rbind") %:%
  foreach(sim_runs = 1:simrunstotal, .combine = "rbind") %dopar%
  {
    seed <- round(runif(1, min = 0, max = 100000)*1000, digits = 0)
    correlation <- simModels$correlation[jj]
    temp <- foreach(distn = 1:nrow(disttable), .packages = c("lavaan", "foreach", "dplyr"), .combine = "rbind") %do%
      {
        data <- lavaan::simulateData(model = simModels$model[jj],
                                     sample.nobs = n, # Number of observations.
                                     skewness = disttable$skewness[[distn]],
                                     kurtosis = disttable$kurtosis[[distn]],
                                     seed = seed, # Set random seed.
                                     empirical = FALSE, # Logical. If TRUE, the implied moments (Mu and Sigma) specify the empirical not population mean and covariance matrix.
                                     return.type = "data.frame"
        )
        
        simuresults <- wrapper_alt(data = data, model = model_est, latent1 = "xi_1", latent2 = "xi_2")
        data.frame (correlation = correlation, 
                    n = n,
                    sim_runs, 
                    seed = seed,
                    type = simModels$type[jj],
                    HTMT_cov = simuresults$htmt_cov,
                    HTMT_cor = simuresults$htmt_cor,
                    HTMT_2_cor = simuresults$htmt_2_cor,
                    MGA = simuresults$mga
                    )
      }
    temp
  }
closeAllConnections()
write.csv2(x = simresults, file = "smallsim_congeneric_mit_seed.csv")

