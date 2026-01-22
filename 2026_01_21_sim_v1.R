.libPaths(c("/home/jab49wd/R-projects/R/library", .libPaths()))

library(lavaan)
library(semTools)
library(doParallel)
library(foreach)
library(dplyr)

source("setupcalc.R")

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
    
    out$model_index <- jj
    out$type <- simModels$type[jj]
    out$correlation <- simModels$correlation[jj]
    out$sample_size <- n
    out$seed <- all_seeds[ii]
    
    
    ii <- ii + 1 
    out
  }

closeAllConnections()
write.csv2(simresults, "simresults_calc2.csv")