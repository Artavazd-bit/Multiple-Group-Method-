.libPaths(c("/home/jab49wd/R-projects/R/library", .libPaths()))

library(lavaan)
#library(MASS)
library(semTools)
#library(gdata)
#library(boot)
#library(dplyr)
#library(stats)
library(doParallel)
library(foreach)
library(dplyr)

source("setupsmall.R")


nkernel <- 8
nobs <- c(50, 100)
simrunstotal <- 10

cl <- parallel::makeCluster(nkernel, outfile = "errorcluster.txt")
doParallel::registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("/home/jab49wd/R-projects/R/library"))

simresults <- foreach(jj = 4:5, .packages = c("lavaan", "foreach", "dplyr", "semTools"), .combine = "rbind") %:%
  foreach(n = nobs, .combine = "rbind") %:%
  foreach(sim_runs = 1:simrunstotal, .combine = "rbind") %do% 
  {
    seed  <- round(runif(1, min = 0, max = 100000)*1000, digits = 0)
    correlation <- simModels$correlation[jj]
    data <-  lavaan::simulateData(model = simModels$model[jj],
                                  sample.nobs = n, # Number of observations.
                                  skewness = NULL,
                                  kurtosis = NULL,
                                  seed = seed, # Set random seed.
                                  empirical = FALSE, # Logical. If TRUE, the implied moments (Mu and Sigma) specify the empirical not population mean and covariance matrix.
                                  return.type = "data.frame"
    )
    out <- bs_p_value(data, 
                      nboot = 500, 
                      model = model_est, 
                      latent1 = "xi_1", 
                      latent2 = "xi_2", 
                      model_unconstrained, 
                      model_constrained
                      )
    mga <- data.frame(Model = jj,
               type = simModels$type[jj],
               correlation = correlation, 
               n = n, 
               seed = seed, 
               type = "mga",
               test_stat = out$test_stat$mga,
               p_value = out$p_values$p_mga, 
               decision = out$p_values$p_mga < 0.05
               ) 
    mga$boot <- list(unlist(out$boot["mga", ]))
    print(mga)
    htmt_cov <- data.frame(Model = jj,
                           type = simModels$type[jj],
                           correlation = correlation, 
                           n = n, 
                           seed = seed, 
                           type = "htmt_cov",
                           test_stat = out$test_stat$htmt_cov,
                           p_value = out$p_values$p_htmt_cov, 
                           decision = out$p_values$p_htmt_cov < 0.05
    )
    htmt_cov$boot <- list(unlist(out$boot["htmt_cov", ]))
    print(htmt_cov)
    htmt_cor <-  data.frame(Model = jj,
                            type = simModels$type[jj],
                            correlation = correlation, 
                            n = n, 
                            seed = seed, 
                            type = "htmt_cor", 
                            test_stat = out$test_stat$htmt_cor,
                            p_value = out$p_values$p_htmt_cor,
                            decision = out$p_values$p_htmt_cor < 0.05
    )
    htmt_cor$boot <- list(unlist(out$boot["htmt_cor", ]))
    print(htmt_cor)
    htmt_2 <-  data.frame(Model = jj,
                            type = simModels$type[jj],
                            correlation = correlation, 
                            n = n, 
                            seed = seed, 
                            type = "htmt_2",
                            test_stat = out$test_stat$htmt_2,
                            p_value = out$p_values$p_htmt_2, 
                            decision = out$p_values$p_htmt_2 < 0.05
    )
    htmt_2$boot <- list(unlist(out$boot["htmt_2_cor", ]))
    print(htmt_2)
    con_phi <-  data.frame(Model = jj,
                          type = simModels$type[jj],
                          correlation = correlation, 
                          n = n, 
                          seed = seed, 
                          type = "con_phi",
                          test_stat = out$test_stat$con_phi,
                          p_value = out$p_values$p_con_phi,
                          decision = out$p_values$p_con_phi < 0.05,
                          boot = NA
    )
    print(con_phi)
    avesv <- list(unname(out$FL_avesv))
    Fornelllarcker <-  data.frame(Model = jj,
                           type = simModels$type[jj],
                           correlation = correlation, 
                           n = n, 
                           seed = seed, 
                           type = "FL",
                           p_value = NA,
                           decision = out$FL_dec,
                           boot = NA
    )
    Fornelllarcker$test_stat <- avesv
    print(Fornelllarcker)
    out2 <- rbind(mga, htmt_cov, htmt_cor, htmt_2, con_phi, Fornelllarcker)
    out2
  }

closeAllConnections()
write.csv2(simresults, "simresults_1.csv")
