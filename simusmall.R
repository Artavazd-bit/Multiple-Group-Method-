library(lavaan)
library(doParallel)
library(foreach)
library(dplyr)
library(matlib)

source("setupsmall.R")

nkernel <- 16
nobs <- c(50, 100, 200)
simrunstotal <- 1000
bootruns <- 500
alphavec <- c(0.01, 0.05, 0.1)
# normal: no skewness, no kurtosis, nonnormal: specified skewness and kurtosis 
# values 
disttable <- data.frame(name = c("normal"))
disttable$skewness <- list(NULL)
disttable$kurtosis <- list(NULL)
######################### Monte Carlo Simulation ###############################
# initialising cluster 
cl <- parallel::makeCluster(nkernel, outfile = "errorcluster.txt")
doParallel::registerDoParallel(cl)

simresults <- foreach(jj = 1:nrow(simModels), .packages = c("lavaan", "foreach", "dplyr"), .combine = "rbind", .export = c("jacknife")) %:%
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
        res <- foreach(type = c("delta", "boot", "bcaboot"), .combine = "rbind") %:%
          foreach(indexgamma = 1:length(alphavec), .combine = "rbind") %do%
          {
            HTMT <- simuresults$delta$HTMT
            seHTMT <- simuresults[[type]]$se
            lowerbound <- unname(simuresults[[type]]$lowerbound[indexgamma])
            upperbound <- unname(simuresults[[type]]$upperbound[indexgamma])
            missing <- unname(simuresults[[type]]$missing)
            time <- simuresults[[type]]$time
            row <- data.frame(correlation = correlation
                              , n = n
                              , HTMT = HTMT
                              , seHTMT = seHTMT
                              , sim_runs
                              , datatype = disttable$name[distn]
                              , method = type
                              , alpha = alphavec[indexgamma]
                              , lowerbound = lowerbound
                              , upperbound =  upperbound
                              , coveragecorr = lowerbound < correlation & correlation <  upperbound
                              , coverageone =  lowerbound < 1 & 1 < upperbound
                              , time = time
                              , seed = seed
                              , missing = missing
            )
            row
          }
        res
      }
    temp
  }
closeAllConnections()
