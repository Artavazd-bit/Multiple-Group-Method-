library(lavaan)
library(semTools)
library(doParallel)
library(foreach)
library(dplyr)

source("2026_01_20_setup_v1.R")

model <- "xi_1 =~ 0.7 * x11 + 0.8 * x12 + 0.9 * x13
          xi_2 =~ 0.8 * x21 + 0.6 * x22 + 0.7 * x23
          xi_1 ~~ 1 * xi_1 + 0.5 * xi_2
          
          x11 ~~ 1 * x11 + 0 * x12 + 0 * x13 + 0 * x21 + 0 * x22 + 0 * x23
          x12 ~~ 1 * x12 + 0 * x13 + 0 * x21 + 0 * x22 + 0 * x23
          x13 ~~ 1 * x13 + 0 * x21 + 0 * x22 + 0 * x23
          
          x21 ~~ 1 * x21 + 0 * x22 + 0 * x23
          x22 ~~ 1 * x22 + 0 * x23
          x23 ~~ 1 * x23"

n <- 100
data <- lavaan::simulateData(model = model,
                              sample.nobs = n, # Number of observations.
                              skewness = NULL,
                              kurtosis = NULL,
                              seed = NULL, # Set random seed.
                              empirical = FALSE, # Logical. If TRUE, the implied moments (Mu and Sigma) specify the empirical not population mean and covariance matrix.
                              return.type = "data.frame"
)

tetrad <- function(data, model, latent1 = NULL, latent2 = NULL, scale = FALSE){
  model_df <- lavaanify(model)
  
  if(is.null(latent1) & is.null(latent2) ){
    latentvars <- unique(model_df$lhs[model_df$op == "=~"])
    listind1 <- list(model_df$rhs[model_df$lhs == latentvars[1] & model_df$op == "=~"])
    listind2 <- list(model_df$rhs[model_df$lhs == latentvars[2] & model_df$op == "=~"])
    if(length(latentvars) != 2){
      error("please select two latent variables")
    }
  }else{
    listind1 <- list(model_df$rhs[model_df$lhs == latent1 & model_df$op == "=~"])
    listind2 <- list(model_df$rhs[model_df$lhs == latent2 & model_df$op == "=~"])
  }
  all_indicators <- unlist(list(listind1, listind2)) 
  
  subset_data <- data[, all_indicators]
  if(scale == FALSE){
    cor_subset_data <- cov(subset_data)
  } else { 
    cor_subset_data <- cor(subset_data)
  }
  
  ind <- which( lower.tri(cor_subset_data,diag=F) , arr.ind = TRUE )
  cor_values <- data.frame( col = dimnames(cor_subset_data)[[2]][ind[,2]] ,
                            row = dimnames(cor_subset_data)[[1]][ind[,1]] ,
                            val = cor_subset_data[ ind ] )
  
  cor_values$type[cor_values$col %in% unlist(listind1) & cor_values$row 
                  %in% unlist(listind1)] <- "mono1"
  cor_values$type[cor_values$col %in% unlist(listind2) & cor_values$row 
                  %in% unlist(listind2)] <- "mono2"
  cor_values$type[cor_values$col %in% unlist(listind1) & cor_values$row 
                  %in% unlist(listind2)] <- "het"
  
  mono <- cor_values[cor_values$type %in% c("mono1", "mono2"), ]
  cor_sym <- cor_values %>%
    bind_rows(
      cor_values %>% rename(col = row, row = col)
    )
  
  out <- mono %>% mutate(id = row_number()) %>% 
    cross_join(., ., suffix = c("_1", "_2")) %>% 
    filter(id_1 < id_2) %>% 
    filter(type_1 != type_2) %>% 
    left_join(., cor_sym, by = c("row_1" = "row", "col_2" = "col")) %>% 
    left_join(., cor_sym, by = c("col_1" = "col", "row_2" = "row"))
  
  out$tetrad <- out$val_1 * out$val_2 - out$val.x*out$val.y
  
  out <- out %>% select(- c("type.x", "type.y", "id_1", "id_2"))
  
  return(out$tetrad)
}
  
t <- tetrad(data, model = model)

bootstrap <- function(data, statisticfun, ...,  alpha = 0.05, nboot)
{
  boot <- sapply(1:nboot, function(x) statisticfun(data = dplyr::sample_n(data, nrow(data), replace = TRUE), ...))
  valid_boot <- boot[!is.na(boot)]
  lowerbound <- unname(quantile(valid_boot, probs = alpha/2))
  upperbound <- unname(quantile(valid_boot, probs = 1 - (alpha/2)))
  bootmean <- mean(valid_boot)
  bootsd <- sd(valid_boot)
  return(list(
    boot = valid_boot,
    lowerbound = lowerbound, 
    upperbound = upperbound,
    mean = bootmean, 
    se = bootsd,
    missing = sum(is.na(boot)),
    alpha = alpha
  ))
}


boot <- bootstrap(data, statisticfun = tetrad, model = model, alpha = 0.05, nboot = 500)

tetrads_boot <- sapply(0:8, function(x) { boot$boot[seq_along(boot$boot) %% 9 == x]})

tetrads_lower <- as.data.frame(tetrads_boot) %>% summarise(across(everything(), ~ quantile(.x, probs = 0.05/2)))

all(tetrads_lower < 0) 

tetrad_test <- function(data, model, latent1 = NULL, latent2 = NULL, scale = FALSE, nboot = 500, alpha = 0.05){
  
  tetrads <- tetrad(data, model = model, latent1 = latent1, latent2 = latent2, scale = scale)
  boot <- bootstrap(data, 
                    statisticfun = tetrad, 
                    model = model, 
                    latent1 = latent1, 
                    latent2 = latent2, 
                    scale = scale, 
                    alpha = 0.05, 
                    nboot = nboot)
  tetrads_boot <- sapply(0:(length(tetrads)-1), function(x) { boot$boot[seq_along(boot$boot) %% length(tetrads) == x]})
  tetrads_lower <- as.data.frame(tetrads_boot) %>% summarise(across(everything(), ~ quantile(.x, probs = alpha/2)))
  tetrads_upper <- as.data.frame(tetrads_boot) %>% summarise(across(everything(), ~ quantile(.x, probs = 1 - alpha/2)))
  
  test <- all(tetrads_lower < 0) & all(tetrads_upper > 0) 
  return(list(test, tetrads_lower = tetrads_lower, tetrads_upper = tetrads_upper))
}

test <- tetrad_test(data, model, alpha = 0.05)


