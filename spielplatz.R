

pop_data <-  lavaan::simulateData(model = simModels$model[4],
                              sample.nobs = 50, # Number of observations.
                              skewness = NULL,
                              kurtosis = NULL,
                              seed = 12345, # Set random seed.
                              empirical = TRUE, # Logical. If TRUE, the implied moments (Mu and Sigma) specify the empirical not population mean and covariance matrix.
                              return.type = "data.frame"
)

multigroup(data = pop_data, model = model_unconstrained, latent1 = "xi_1", latent2 = "xi_2")

calchtmt(data = pop_data, model = model_unconstrained, latent1 = "xi_1", latent2 = "xi_2", scale = TRUE, htmt2 = TRUE)
