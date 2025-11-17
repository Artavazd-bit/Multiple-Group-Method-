library(lavaan)
library(doParallel)
library(foreach)
library(dplyr)
library(matlib)



coefs <- 1
corr <- c(0.7, 0.8, 0.9, 0.95, 1)
param <- expand.grid(loading1 = coefs, loading2 = coefs, correlation = corr)
simModels <- foreach(i = 1:nrow(param), .combine = "rbind") %do%
  {
    simCommonFactor <- 
      paste(
        paste("xi_1 =~ ",param$loading1[i], "*x11 + ",param$loading1[i],"*x12 + ",param$loading1[i], "*x13"),"\n"
        , paste("xi_2 =~ ",param$loading2[i], "*x21 + ",param$loading2[i],"*x22 + ",param$loading2[i], "*x23"), "\n"
        , paste("xi_1 ~~ 1*xi_1 + ", param$correlation[i], "*xi_2"),"\n"
        , "xi_2 ~~ 1*xi_2 \n"
        , paste("x11 ~~", 0.6, "*x11 + 0*x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x12 ~~", 0.5, "*x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x13 ~~", 0.2, "*x13 + 0*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x21 ~~", 0.6, "*x21 + 0*x22 + 0*x23"),"\n"
        , paste("x22 ~~", 0.5, "*x22 + 0*x23"),"\n"
        , paste("x23 ~~", 0.2, "*x23"), "\n"
        , paste("x11 ~ 0*1"), "\n"
        , paste("x12 ~ 0*1"), "\n"
        , paste("x13 ~ 0*1"), "\n"
        , paste("x21 ~ 0*1"), "\n"
        , paste("x22 ~ 0*1"), "\n"
        , paste("x23 ~ 0*1"), "\n"
      )
    save <- data.frame(
      loading_1 = param$loading1[i],
      loading_2 =  param$loading2[i],
      correlation = param$correlation[i],
      model = simCommonFactor
    )
    save
    #rm(save, simCommonFactor, i)
  }
simModels

rm(coefs, corr, param, simCommonFactor, save, i)

model <- "xi_1 =~  0.8 *x11 +  0.9 *x12 +  0.7 *x13 
xi_2 =~  0.6 *x21 +  0.9 *x22 +  0.8 *x23 
xi_1 ~~ 1*xi_1 +  0.8 *xi_2 
xi_2 ~~ 1*xi_2 
x11 ~~ 0.6 *x11 + 0*x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23 
x12 ~~ 0.5 *x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23 
x13 ~~ 0.2 *x13 + 0*x21 + 0*x22 + 0*x23 
x21 ~~ 0.6 *x21 + 0*x22 + 0*x23 
x22 ~~ 0.5 *x22 + 0*x23 
x23 ~~ 0.2 *x23 
x11 ~ 0*1  
x12 ~ 0*1 
x13 ~ 0*1 
x21 ~ 0*1 
x22 ~ 0*1 
x23 ~ 0*1 "


model2 <- "xi_1 =~  0.8 *x11 +  0.9 *x12 +  0.7 *x13 
xi_2 =~  0.6 *x21 +  0.9 *x22 +  0.8 *x23 
xi_1 ~~ 1*xi_1 +  0.8 *xi_2 
xi_2 ~~ 1*xi_2 
x11 ~~ (1-0.8^2) *x11 + 0*x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23 
x12 ~~ (1-0.9^2) *x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23 
x13 ~~ (1-0.7^2) *x13 + 0*x21 + 0*x22 + 0*x23 
x21 ~~ (1-0.6^2) *x21 + 0*x22 + 0*x23 
x22 ~~ (1-0.9^2) *x22 + 0*x23 
x23 ~~ (1-0.8^2) *x23 
x11 ~ 0*1   
x12 ~ 0*1 
x13 ~ 0*1 
x21 ~ 0*1 
x22 ~ 0*1 
x23 ~ 0*1 "
################################################################################
model_est<- '
              #  latent variables
                xi_1 =~ NA*x11 + x12 + x13
                xi_2 =~ NA*x21 + x22 + x23 
                
                xi_1 ~~ xi_2
                xi_1 ~~ 1*xi_1
                xi_2 ~~ 1*xi_2
            ' 
################################################################################

set_seed <- 12345
n <- 10000

# simModels$model[2] gibt eine Korrelation gleich 0.8
data <- lavaan::simulateData(model = model,
                             sample.nobs = n, # Number of observations.
                             skewness = NULL,
                             kurtosis = NULL,
                             seed = set_seed, # Set random seed.
                             empirical = TRUE, # if TRUE dann sind empirical gleich pop values
                             return.type = "data.frame"
)
################################################################################
out <- cfa(model = model_est, data = data)
summary(out, standardized = TRUE)
################################################################################
# communalities passen noh nicht, weiter ausprobieren
S <- cor(data)
Sstar <- S
dimcount <- 1:dim(Sstar)[1]
communialities <- sapply(dimcount, function(x) {sum(combn(Sstar[x,-x], m = 2, FUN = prod))/sum(Sstar[-x,-x][upper.tri(Sstar[-x,-x])])}) 

S1 <- cor(data[,1:3])
Sstar1 <- S1
dimcount1 <- 1:dim(Sstar1)[1]
communialities1 <- sapply(dimcount1, function(x) {sum(combn(Sstar1[x,-x], m = 2, FUN = prod))/sum(Sstar1[-x,-x][upper.tri(Sstar1[-x,-x])])}) 

S2 <- cor(data[,4:6])
Sstar2 <- S2
dimcount2 <- 1:dim(Sstar2)[1]
communialities2 <- sapply(dimcount2, function(x) {sum(combn(Sstar2[x,-x], m = 2, FUN = prod))/sum(Sstar2[-x,-x][upper.tri(Sstar2[-x,-x])])}) 

diag(Sstar) <- c(communialities1, communialities2)
# lavInspect(out, "communalities")

model_df <- lavaanify(model_est)
latent1 <- "xi_1"
latent2 <- "xi_2"
listind1 <- list(model_df$rhs[model_df$lhs == latent1 & model_df$op == "=~"])
listind2 <- list(model_df$rhs[model_df$lhs == latent2 & model_df$op == "=~"])
pattern_matrix <- S[,c(1,2)] * 0 
colnames(pattern_matrix) <- c(latent1, latent2)
pattern_matrix[unlist(listind1),latent1] <- 1 
pattern_matrix[unlist(listind2),latent2] <- 1 

pattern_matrix <- t(pattern_matrix)

Fzero <- pattern_matrix %*% Sstar 
Lzero <- Fzero %*% t(pattern_matrix)

D <- diag(diag(Lzero))
Dinvsqrt <- diag(1/sqrt(diag(D)))

L <- Dinvsqrt %*% Lzero %*% Dinvsqrt

Fmatrix <- t(Fzero) %*% Dinvsqrt  
################################################################################
source("setupmgm.R")


multigroup(data = data, model = model_est, latent1 = "xi_1", latent = "xi_2")
calchtmt_alt(data = data, model = model_est, latent1 = "xi_1", latent = "xi_2", scale = FALSE, htmt2 = FALSE)

