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
    
    out <- tryCatch(
      {
        wrapper_function(data = data, model = model_est, latent1 = "xi_1", latent2 = "xi_2", 
                         model_constrained = model_constrained, 
                         model_unconstrained = model_unconstrained)
      },
      error = function(e) {
        cat("Error in iteration jj=", jj, " sim_runs=", sim_runs, " seed=",  all_seeds[ii], "\n")
        cat("Error message:", conditionMessage(e), "\n")
        # Return a valid 6-row data frame with NAs
        data.frame(
          test = c("MGA", "HTMT_COV", "HTMT_COR", "HTMT_2", "CON_PHI", "FL"),
          stat = NA, constrained = NA, unconstrained = NA,
          FL_1 = NA, FL_2 = NA, warning = NA, error = conditionMessage(e)
        )
      }
    )
    out$type <- simModels$type[jj]
    out$correlation <- simModels$correlation[jj]
    out$seed <- all_seeds[ii]
    ii <- ii + 1 
    out
  }

closeAllConnections()
write.csv2(simresults, "simresults_calc.csv")