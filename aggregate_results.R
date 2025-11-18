library(lavaan)
library(doParallel)
library(foreach)
library(dplyr)
library(matlib)
simresults2 <- simresults

simresults2$mgamissing <- as.numeric(is.na(simresults2$MGA))
simresults2$htmt_cov_missing <- as.numeric(is.na(simresults2$HTMT_cov))
simresults2$htmt_cor_missing <- as.numeric(is.na(simresults2$HTMT_cor))
simresults2$htmt_2_cor_missing <- as.numeric(is.na(simresults2$HTMT_2_cor))


missingag <- simresults2 %>% group_by(correlation, n, type) %>%
  summarize(mga_missing = mean(mgamissing),
            htmt_cov_missing = mean(htmt_cov_missing),
            htmt_cor_missing = mean(htmt_cor_missing),
            htmt_2_cor_missing = mean(htmt_2_cor_missing)
  )

resultsag <- simresults %>% 
  group_by(correlation, n, loading1, loading2) %>%
  summarize(htmt_cov_mean = mean(HTMT_cov),
            htmt_cor_mean = mean(HTMT_cor),
            htmt_2_cor_mean = mean(HTMT_2_cor),
            mga_mean = mean(MGA)
  )


missingag_sort <- missingag %>% arrange(desc(mga_missing), desc(htmt_cov_missing), desc(htmt_cor_missing), desc(htmt_2_cor_missing))

######### find missing #############
seed_missing_mga <- simresults[simresults$correlation == 0 & simresults$n == 25 & simresults$type == "harsh" & is.na(simresults$MGA),]

######### missing #################
data <- lavaan::simulateData(model = simModels$model[1],
                             sample.nobs = 25, # Number of observations.
                             skewness = NULL,
                             kurtosis = NULL,
                             seed = 54564242, # Set random seed.
                             empirical = FALSE, # if TRUE dann sind empirical gleich pop values
                             return.type = "data.frame"
)

multigroup(data = data, mode = model_est, latent1 = "xi_1", latent2 = "xi_2")






