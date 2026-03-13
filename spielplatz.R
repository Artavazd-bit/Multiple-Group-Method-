
source("2026_01_20_setup_v1.R")
install.packages("Matrix")
library(Matrix)


pop_data <-  lavaan::simulateData(model = simModels$model[5],
                              sample.nobs = 10000, # Number of observations.
                              skewness = NULL,
                              kurtosis = NULL,
                              seed = 12345, # Set random seed.
                              empirical = FALSE, # Logical. If TRUE, the implied moments (Mu and Sigma) specify the empirical not population mean and covariance matrix.
                              return.type = "data.frame"
)

multigroup(data = pop_data, model = model_unconstrained, latent1 = "xi_1", latent2 = "xi_2")

calchtmt(data = pop_data, model = model_unconstrained, latent1 = "xi_1", latent2 = "xi_2", scale = TRUE, htmt2 = TRUE)

model_x11 <- lm(x11 ~ x21 + x22 + x23, pop_data)
rsq_x11 <- summary(model_x11)$r.squared

model_testx11 <- lm(x11 ~ x12 + x13, pop_data)
rsq_x11_own <- summary(model_testx11)$r.squared

rsq_x11 / rsq_x11_own 

model_x11_x21 <- lm(x11 ~ x21, pop_data)
summary(model_x11_x21)

model_x12 <- lm(x12 ~ x21 + x22 + x23, pop_data)
rsq_x12 <- summary(model_x12)$r.squared

model_x13 <- lm(x13 ~ x21 + x22 + x23, pop_data)
rsq_x13 <- summary(model_x13)$r.squared

l1 <- c(rsq_x11, rsq_x12, rsq_x13)


model_x21 <- lm(x21 ~ x11 + x12 + x13, pop_data)
rsq_x21 <- summary(model_x21)$r.squared

model_x22 <- lm(x22 ~ x11 + x12 + x13, pop_data)
rsq_x22 <- summary(model_x22)$r.squared

model_x23 <- lm(x23 ~ x11 + x12 + x13, pop_data)
rsq_x23 <- summary(model_x23)$r.squared

l2 <- c(rsq_x21, rsq_x22, rsq_x23)
mean(l1)
mean(l2)
maxl1 <- max(l1)
maxl2 <- max(l2)
1/(1- mean(l1))
1/(1- exp(mean(log(l1))))
1/(1- maxl1)
1/(1- mean(l2))
1/(1- exp(mean(log(l2))))
1/(1-maxl2)


colnames(pop_data)

model_x11x21 <- lm(x11 ~ x21, pop_data)
rsqx11x21 <- summary(model_x11x21)$r.squared

model_x11x22 <- lm(x11 ~ x22, pop_data)
rsqx11x22 <- summary(model_x11x22)$r.squared

model_x11x23 <- lm(x11 ~ x23, pop_data)
rsqx11x23 <- summary(model_x11x23)$r.squared

x11 <- c(rsqx11x21, rsqx11x22, rsqx11x23)


model_x12x21 <- lm(x12 ~ x21, pop_data)
rsqx12x21 <- summary(model_x12x21)$r.squared

model_x12x22 <- lm(x12 ~ x22, pop_data)
rsqx12x22 <- summary(model_x12x22)$r.squared

model_x12x23 <- lm(x12 ~ x23, pop_data)
rsqx12x23 <- summary(model_x12x23)$r.squared

x12 <- c(rsqx12x21, rsqx12x22, rsqx12x23)

lm(c(x11, x12, x13) ~ x21 + x22 + x23)



X <- as.matrix(pop_data)
P <- matrix(c(1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1), nrow = 2, ncol = 6)


test <- P %*% (t(X) %*% X) %*% t(P)


angle <- function(a, b) {
  cos_theta <- sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
  acos(cos_theta)  # Ergebnis in Radiant
}

area <- function(a, b) {
  abs(a[1] * b[2] - a[2] * b[1])
}

angle(test[,1], test[,2])
area(test[,1], test[,2])


model_df <- lavaanify(model)
listind1 <- list(model_df$rhs[model_df$lhs == latent1 & model_df$op == "=~"])
listind2 <- list(model_df$rhs[model_df$lhs == latent2 & model_df$op == "=~"])

S <- cor(data[,c(unlist(listind1), unlist(listind2))])
Sstar <- S

corlist <- list(cor(data[,unlist(listind1)]), cor(data[,unlist(listind2)]) ) 
communalities <- lapply(corlist, function(corlist) {
  sapply(1:nrow(corlist), function(x) {sum(combn(corlist[x,-x], m = 2, FUN = prod))/sum(corlist[-x,-x][upper.tri(corlist[-x,-x])])})
})

diag(Sstar) <- unlist(communalities)

Sstar
rankMatrix(Sstar)


model_est <- "xi_1 =~ 0.7 * x11 + 0.8 * x12 + 0.9 * x13
          xi_2 =~ 0.8 * x21 + 0.6 * x22 + 0.7 * x23
          xi_1 ~~ 1 * xi_1 + 1 * xi_2
          
          x11 ~~ 1 * x11 + 0 * x12 + 0 * x13 + 0 * x21 + 0 * x22 + 0 * x23
          x12 ~~ 1 * x12 + 0 * x13 + 0 * x21 + 0 * x22 + 0 * x23
          x13 ~~ 1 * x13 + 0 * x21 + 0 * x22 + 0 * x23
          
          x21 ~~ 1 * x21 + 0 * x22 + 0 * x23
          x22 ~~ 1 * x22 + 0 * x23
          x23 ~~ 1 * x23"


data <-  lavaan::simulateData(model = model_est,
                                  sample.nobs = 10000, # Number of observations.
                                  skewness = NULL,
                                  kurtosis = NULL,
                                  seed = 12345, # Set random seed.
                                  empirical = TRUE, # Logical. If TRUE, the implied moments (Mu and Sigma) specify the empirical not population mean and covariance matrix.
                                  return.type = "data.frame"
)

cor(data)
cov(data) * 9999/10000

data2 <- data * sqrt(9999/10000)

cov(data1)
cov(data2)

model_df <- lavaanify(model_unconstrained)


compute_2x2_minors <- function(M) {
  n <- nrow(M)
  rows <- combn(n, 2)
  cols <- combn(n, 2)
  
  result <- matrix(NA, nrow = ncol(rows), ncol = ncol(cols))
  rownames(result) <- apply(rows, 2, \(r) paste0("r", r[1], r[2]))
  colnames(result) <- apply(cols, 2, \(c) paste0("c", c[1], c[2]))
  
  for (i in seq_len(ncol(rows))) {
    for (j in seq_len(ncol(cols))) {
      sub <- M[rows[, i], cols[, j]]
      result[i, j] <- sub[1, 1] * sub[2, 2] - sub[1, 2] * sub[2, 1]
    }
  }
  
  result
}


test <- compute_2x2_minors(cov(data2))

m <- cov(data2)

tetrad1 <- m["x11","x12"] * m["x21", "x22"] - m["x11", "x21"] * m["x12", "x22"]

tetrad

n <- nrow(m)
rows <- combn(n, 2)

