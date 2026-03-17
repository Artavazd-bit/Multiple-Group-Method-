library(lavaan)
library(semTools)
library(doParallel)
library(foreach)
library(dplyr)
library(numDeriv)

source("Simu2/2026_03_07_setup_v2.R")

cov <- cov(data_1)
x <- cov[lower.tri(cov)]

matrix <- matrix(nrow = 6, ncol = 6)

matrix[lower.tri(matrix)] <- x

matrix[upper.tri(matrix)] <- t(matrix)[upper.tri(matrix)]

ind <- which( lower.tri(matrix,diag=F) , arr.ind = TRUE )

cor_values <- data.frame( col = dimnames(matrix)[[2]][ind[,2]] ,
                          row = dimnames(matrix)[[1]][ind[,1]] ,
                          val = matrix[ ind ] )

test <- foreach(i = 1:1000, .combine = "rbind") %do%{
  data_1 <- lavaan::simulateData(model = model_cor1,
                                 sample.nobs = n, # Number of observations.
                                 skewness = NULL,
                                 kurtosis = NULL,
                                 seed = NULL, # Set random seed.
                                 empirical = FALSE, # Logical. If TRUE, the implied moments (Mu and Sigma) specify the empirical not population mean and covariance matrix.
                                 return.type = "data.frame"
  )
  t <- tetrad(data_1, model_unconstrained)
  t$tetrad
}

testcov <- cov(test)



model <- "xi_1 =~ 0.7 * x11 + 0.8 * x12 
          xi_2 =~ 0.8 * x21 + 0.6 * x22 
          xi_1 ~~ 1 * xi_1 + 1 * xi_2
          
          x11 ~~ 1 * x11 + 0 * x12  + 0 * x21 + 0 * x22 
          x12 ~~ 1 * x12  + 0 * x21 + 0 * x22 
          
          x21 ~~ 1 * x21 + 0 * x22 
          x22 ~~ 1 * x22 "
