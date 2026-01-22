.libPaths(c("/home/jab49wd/R-projects/R/library", .libPaths()))

library(lavaan)
library(semTools)
library(doParallel)
library(foreach)
library(dplyr)

source("2026_01_20_setup_v1.R")

nkernel <- 32
nobs <- c(50, 100, 200, 1000)
simrunstotal <- 10000
all_seeds <- sample(1:1e9, size = nrow(simModels) * length(nobs) * simrunstotal)
ii <- 1
cl <- parallel::makeCluster(nkernel, outfile = "errorcluster.txt")
doParallel::registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("/home/jab49wd/R-projects/R/library"))

simresults <- foreach(jj = 1:nrow(simModels), .packages = c("lavaan", "foreach", "dplyr", "semTools"), .combine = "rbind") %:%
  foreach(n = nobs, .combine = "rbind") %:%
  foreach(sim_runs = 1:simrunstotal, .combine = "rbind") %do% 
  {
    correlation <- simModels$correlation[jj]
    data <-  lavaan::simulateData(model = simModels$model[jj],
                                  sample.nobs = n, # Number of observations.
                                  skewness = NULL,
                                  kurtosis = NULL,
                                  seed = all_seeds[ii], # Set random seed.
                                  empirical = FALSE, # Logical. If TRUE, the implied moments (Mu and Sigma) specify the empirical not population mean and covariance matrix.
                                  return.type = "data.frame"
    )
    out <- run_methods(data = data,
                       model_unconstrained = model_unconstrained, 
                       latent1 = "xi_1", 
                       latent2 = "xi_2", 
                       model_constrained = model_constrained)
    
    
    out$model_index <- jj
    out$type <- simModels$type[jj]
    out$correlation <- simModels$correlation[jj]
    out$sample_size <- n
    out$seed <- all_seeds[ii]
    
    
    ii <- ii + 1 
    out
  }

closeAllConnections()
write.csv2(simresults, "simres_v1.csv")