library(lavaan)
library(semTools)
library(doParallel)
library(foreach)
library(dplyr)
library(numDeriv)

source("Simu2/2026_03_07_setup_v2.R")

model_cor1 <- simModels$model[4]
model_cor03 <- simModels$model[2]

n <- 10000
data_1 <- lavaan::simulateData(model = model_cor1,
                                       sample.nobs = n, # Number of observations.
                                       skewness = NULL,
                                       kurtosis = NULL,
                                       seed = NULL, # Set random seed.
                                       empirical = TRUE, # Logical. If TRUE, the implied moments (Mu and Sigma) specify the empirical not population mean and covariance matrix.
                                       return.type = "data.frame"
)
data_03 <- lavaan::simulateData(model = model_cor03,
                                       sample.nobs = n, # Number of observations.
                                       skewness = NULL,
                                       kurtosis = NULL,
                                       seed = NULL, # Set random seed.
                                       empirical = TRUE, # Logical. If TRUE, the implied moments (Mu and Sigma) specify the empirical not population mean and covariance matrix.
                                       return.type = "data.frame"
)



t1 <- tetrad(data_1, model_unconstrained)
t3 <- tetrad(data_03, model_unconstrained)

tetrad_test(data_1, model = model_unconstrained)
tetrad_test(data_03, model = model_unconstrained)

tetrad_outer(data_1, model_unconstrained)


out3 <- tetrad_outer(data_03, model_unconstrained)

covcov3 <- calcovcov(data_03)

covtetrads3 <- out3$grad_tetrads %*% (covcov3 / nrow(data_03) )  %*% t(out3$grad_tetrads)
T <- nrow(data_03)* t(out3$tetrads) %*% solve(covtetrads3) %*% out3$tetrads 
