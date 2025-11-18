multigroup <- function(data , model, latent1, latent2){
  model_df <- lavaanify(model_est)
  listind1 <- list(model_df$rhs[model_df$lhs == latent1 & model_df$op == "=~"])
  listind2 <- list(model_df$rhs[model_df$lhs == latent2 & model_df$op == "=~"])
  
  S <- cor(data[,c(unlist(listind1), unlist(listind2))])
  Sstar <- S
  
  corlist <- list(cor(data[,unlist(listind1)]), cor(data[,unlist(listind2)]) ) 
  communalities <- lapply(corlist, function(corlist) {
    sapply(1:nrow(corlist), function(x) {sum(combn(corlist[x,-x], m = 2, FUN = prod))/sum(corlist[-x,-x][upper.tri(corlist[-x,-x])])})
  })
  
  diag(Sstar) <- unlist(communalities)
  
  pattern_matrix <- S[,c(1,2)] * 0 
  colnames(pattern_matrix) <- c(latent1, latent2)
  pattern_matrix[unlist(listind1),latent1] <- 1 
  pattern_matrix[unlist(listind2),latent2] <- 1 
  pattern_matrix <- t(pattern_matrix)
  
  Fzero <- pattern_matrix %*% Sstar 
  Lzero <- Fzero %*% t(pattern_matrix)
  
  Dinvsqrt <- diag(1/sqrt(diag(Lzero)))
  
  L <- Dinvsqrt %*% Lzero %*% Dinvsqrt
  
  Fmatrix <- t(Fzero) %*% Dinvsqrt 
  
  return(L[2,1])
}
################################################################################
calchtmt <- function(data, model, latent1, latent2, scale, htmt2){
  
  model_df <- lavaanify(model)
  
  listind1 <- list(model_df$rhs[model_df$lhs == latent1 & model_df$op == "=~"])
  listind2 <- list(model_df$rhs[model_df$lhs == latent2 & model_df$op == "=~"])
  
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
  
  K_i <- length(unlist(listind1))
  K_j <- length(unlist(listind2))
  
  tryCatch({
    if(htmt2 == FALSE){
      A = 1/(K_i*K_j) * sum(cor_values$val[cor_values$type == "het"])
      B = 2/(K_i*(K_i-1)) *  sum(cor_values$val[cor_values$type == "mono1"]) 
      C = 2/(K_j*(K_j-1)) *  sum(cor_values$val[cor_values$type == "mono2"]) 
      HTMT <- A / ((B*C)^(1/2))
    }
    else if(htmt2 == TRUE){
      A =  prod(cor_values$val[cor_values$type == "het"])^(1/(K_i*K_j))
      B =  prod(cor_values$val[cor_values$type == "mono1"])^(2/(K_i*(K_i-1))) 
      C =  prod(cor_values$val[cor_values$type == "mono2"])^(2/(K_j*(K_j-1))) 
      HTMT <- A / ((B*C)^(1/2))
    }
    else{
      print("ERROR")
    }
    return(HTMT)
  } , error = function(e){
    cat("An error occurred:", e$message, "\n")
    return(NA)
  }, warning = function(w){
    cat("A warning occurred:", w$message, "\n")
    return(NA)
  })
}
###############################################################################
wrapper_alt <- function(data, model, latent1, latent2){
  htmt_cov <- calchtmt(data = data, model = model, latent1 = latent1, latent2 = latent2, scale = FALSE, htmt2 = FALSE)
  htmt_cor <- calchtmt(data = data, model = model, latent1 = latent1, latent2 = latent2, scale = TRUE, htmt2 = FALSE)
  htmt_2_cor <- calchtmt(data = data, model = model, latent1 = latent1, latent2 = latent2, scale = TRUE, htmt2 = TRUE)
  
  mga <- multigroup(data = data, model = model, latent1 = latent1, latent2 = latent2)
  
  list(mga = mga, htmt_cov = htmt_cov, htmt_cor = htmt_cor, htmt_2_cor = htmt_2_cor)
}

################################################################################
## Defining the models
#coefs <- c(-0.9, -0.5, 0.01, 0.5, 0.9)
#corr <- c(0.7, 1)
#param <- expand.grid(loading1 = coefs, loading2 = coefs, correlation = corr)
#simModels2 <- foreach(i = 1:nrow(param), .combine = "rbind") %do%
#  {
#    simCommonFactor <- 
#      paste(
#        paste("xi_1 =~ ",param$loading1[i], "*x11 + ",param$loading1[i],"*x12 + ",param$loading1[i], "*x13"),"\n"
#        , paste("xi_2 =~ ",param$loading2[i], "*x21 + ",param$loading2[i],"*x22 + ",param$loading2[i], "*x23"), "\n"
#        , paste("xi_1 ~~ 1*xi_1 + ", param$correlation[i], "*xi_2"),"\n"
#        , "xi_2 ~~ 1*xi_2 \n"
#        , paste("x11 ~~", 0.6, "*x11 + 0*x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23"),"\n"
#        , paste("x12 ~~", 0.5, "*x12 + 0*x13 + 0*x21 + 0*x22 + 0*x23"),"\n"
  #       , paste("x13 ~~", 0.2, "*x13 + 0*x21 + 0*x22 + 0*x23"),"\n"
  #       , paste("x21 ~~", 0.6, "*x21 + 0*x22 + 0*x23"),"\n"
  #       , paste("x22 ~~", 0.5, "*x22 + 0*x23"),"\n"
  #       , paste("x23 ~~", 0.2, "*x23"), "\n"
  #       , paste("x11 ~ 0*1"), "\n"
  #       , paste("x12 ~ 0*1"), "\n"
  #       , paste("x13 ~ 0*1"), "\n"
  #       , paste("x21 ~ 0*1"), "\n"
  #       , paste("x22 ~ 0*1"), "\n"
  #       , paste("x23 ~ 0*1"), "\n"
  #     )
  #   save <- data.frame(
  #     loading_1 = param$loading1[i],
  #     loading_2 =  param$loading2[i],
  #     correlation = param$correlation[i],
  #     model = simCommonFactor
  #   )
  #   save
  #   #rm(save, simCommonFactor, i)
  # }
#simModels2

################################################################################
#coefs <- c(-0.9, -0.5, 0.01, 0.5, 0.9)
corr <- c(0, 0.3, 0.7, 1)
param <- expand.grid(correlation = corr)
simModels_harsh <- foreach(i = 1:nrow(param), .combine = "rbind") %do%
  {
    simCommonFactor <- 
      paste(
        paste("xi_1 =~ 0.7*x11 + (-0.5)*x12 + 0.8*x13"),"\n"
        , paste("xi_2 =~ 0.9*x21 + (-0.6)*x22 + (-0.8)*x23"), "\n"
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
      type = "harsh",
      correlation = param$correlation[i],
      model = simCommonFactor
    )
    save
    #rm(save, simCommonFactor, i)
  }
simModels_harsh

corr <- c(0, 0.3, 0.7, 1)
param <- expand.grid(correlation = corr)
simModels_mild <- foreach(i = 1:nrow(param), .combine = "rbind") %do%
  {
    simCommonFactor <- 
      paste(
        paste("xi_1 =~ 0.7*x11 + (0.5)*x12 + 0.8*x13"),"\n"
        , paste("xi_2 =~ 0.9*x21 + (0.6)*x22 + (0.8)*x23"), "\n"
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
      type = "mild",
      correlation = param$correlation[i],
      model = simCommonFactor
    )
    save
    #rm(save, simCommonFactor, i)
  }

simModels <- rbind(simModels_harsh, simModels_mild)
rm(simModels_harsh, simModels_mild, param, save, corr)
################################################################################
model_est<- '
              #  latent variables
                xi_1 =~ x11 + x12 + x13
                xi_2 =~ x21 + x22 + x23 
                
                xi_1 ~~ xi_2
              ' 
