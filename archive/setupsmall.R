multigroup <- function(data , model, latent1, latent2){
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
  
  pattern_matrix <- S[,c(1,2)] * 0 
  colnames(pattern_matrix) <- c(latent1, latent2)
  pattern_matrix[unlist(listind1),latent1] <- 1 
  pattern_matrix[unlist(listind2),latent2] <- 1 
  pattern_matrix <- t(pattern_matrix)
  
  Fzero <- pattern_matrix %*% Sstar 
  Lzero <- Fzero %*% t(pattern_matrix)
  
  warningmes <- NA
  errormes <- NA
  
  tryCatch(
    {
      Dinvsqrt <- diag(1/sqrt(diag(Lzero)))
      L <- Dinvsqrt %*% Lzero %*% Dinvsqrt
      Fmatrix <- t(Fzero) %*% Dinvsqrt
      return(list(out = L[2,1], error = errormes, warning = warningmes)) 
    },warning = function(w)
    {
      warningmes <<- conditionMessage(w)
      return(list(out =  NA, error = errormes, warning = warningmes))
    }, error = function(e)
    {
      errormes <<- conditionMessage(e)
      return(list(out = NA, error = errormes, warning = warningmes))
    }
  )
  
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
  
  warningmes <- NA
  errormes <- NA
  tryCatch({
    if(htmt2 == FALSE){
      A = 1/(K_i*K_j) * sum(cor_values$val[cor_values$type == "het"])
      B = 2/(K_i*(K_i-1)) *  sum(cor_values$val[cor_values$type == "mono1"]) 
      C = 2/(K_j*(K_j-1)) *  sum(cor_values$val[cor_values$type == "mono2"]) 
      HTMT <- A / ((B*C)^(1/2))
    }
    else if(htmt2 == TRUE){
      A =  (prod(cor_values$val[cor_values$type == "het"]))^(1/(K_i*K_j))
      B =  (prod(cor_values$val[cor_values$type == "mono1"]))^(2/(K_i*(K_i-1))) 
      C =  (prod(cor_values$val[cor_values$type == "mono2"]))^(2/(K_j*(K_j-1))) 
      HTMT <- A / ((B*C)^(1/2))
    }
    else{
      print("ERROR")
    }
    if(is.na(A)){
      errormes <- paste(errormes, "Nominator is NaN") 
    }
    if(is.na(B) | is.na(C)){
      errormes <- paste(errormes, "Denominator is NaN")
    }
    if(B*C < 0){
      errormes <- paste(errormes, "Negative Number under root")
    }
    return(list(out = HTMT, error = errormes, warning = NA)) 
  } , error = function(e){
    errormes <<- conditionMessage(e)
    return(list(out = NA, error = errormes, warning = warningmes) ) 
  }, warning = function(w){
    warningmes <<- conditionMessage(w)
    return(list(out = HTMT, error = errormes, warning = warningmes))
  })
}
###############################################################################
# wrapper_alternativ <- function(data, model, latent1, latent2){
#   htmt_cov <- calchtmt(data = data, model = model, latent1 = latent1, latent2 = latent2, scale = FALSE, htmt2 = FALSE)
#   htmt_cor <- calchtmt(data = data, model = model, latent1 = latent1, latent2 = latent2, scale = TRUE, htmt2 = FALSE)
#   htmt_2_cor <- calchtmt(data = data, model = model, latent1 = latent1, latent2 = latent2, scale = TRUE, htmt2 = TRUE)
#   
#   mga <- multigroup(data = data, model = model, latent1 = latent1, latent2 = latent2)
#   
#   list(mga = mga, htmt_cov = htmt_cov, htmt_cor = htmt_cor, htmt_2_cor = htmt_2_cor)
# }


constrained_phi <- function(model_constrained, model_unconstrained, data){

  constrained_result <- tryCatch(
    {
      warn_msg <- NA
      fit <- withCallingHandlers(
        sem(model_constrained, data = data),
        warning = function(w) {
          warn_msg <<- conditionMessage(w)
        }
      )
      list(fit = fit, warning = warn_msg, error = NA)
    },
    error = function(e){
      list(fit = NA, warning = NA, error = conditionMessage(e))
    }
  )
  
  unconstrained_result <- tryCatch(
    {
    warn_msg <- NA
    fit <- withCallingHandlers(
      sem(model_unconstrained, data = data),
      warning = function(w) {
        warn_msg <<- conditionMessage(w)
      }
    )
    list(fit = fit, warning = warn_msg, error = NA)
    },
    error = function(e){
      list(fit = NA, warning = NA, error = conditionMessage(e))
    }
  )
  if (is.na(constrained_result$error) & is.na(unconstrained_result$error)) {
    test_stat <- constrained_result$fit@test$standard$stat - unconstrained_result$fit@test$standard$stat
    p <- 1 - pchisq(test_stat, df = 3)
    constrained_chi <- constrained_result$fit@test$standard$stat
    unconstrained_chi <- unconstrained_result$fit@test$standard$stat
    }else {
    test_stat <- NA
    p <- NA
    constrained_chi <- NA
    unconstrained_chi <- NA
    }
  if(is.null(constrained_chi)){
    constrained_chi <- NA
  }
  if(is.null(unconstrained_chi)){
    unconstrained_chi <- NA
  }
  
  return(list(
    p = p, 
    test_stat = test_stat, 
    constrained = constrained_chi, 
    unconstrained = unconstrained_chi,
    constrained_warning = constrained_result$warning,
    constrained_error = constrained_result$error,
    unconstrained_warning = unconstrained_result$warning,
    unconstrained_error = unconstrained_result$error
  ))
 }

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

bootstrap_p_value <- function(data = data, statisticfun, ..., nboot)
{
  boot <- sapply(1:nboot, function(x) statisticfun(data = dplyr::sample_n(data, nrow(data), replace = TRUE), ...))
  valid_boot <- boot[!is.na(boot)]
  p <- length(valid_boot[valid_boot > 1]) / length(valid_boot)
}

fornell_larcker <- function(data, model, latent1, latent2){
  result <- tryCatch(
    {
      warn_msg <- NA
      fit <- withCallingHandlers(
        sem(model = model, data = data),
        warning = function(w) {
          warn_msg <<- conditionMessage(w)
        }
      )
      list(fit = fit, warning = warn_msg, error = NA)
    },
    error = function(e){
      list(fit = NA, warning = NA, error = conditionMessage(e))
    }
  ) 
  
  if(is.na(result$warning) & is.na(result$error)){
    ave <- AVE(result$fit)
    ave_sel <- ave[c(latent1, latent2)]
    cormatrix <- lavInspect(fit, "cor.lv")
    avesv <- ave_sel / cormatrix[latent1, latent2] 
    decision <- any(avesv < 1)
  }else{
    decision <- NA
    avesv <- NA
  }
  
  return(list(dec = decision, 
              ave_sv = avesv,
              warning = result$warning,
              error = result$error))
}

wrap <- function(data, model_est, latent1, latent2){
  mga_value <- multigroup(data = data, model = model_est, 
                          latent1 = latent1, latent2 = latent2)
  htmt_cov <- calchtmt(data = data, model = model_est, latent1 = latent1, latent2 = latent2, scale = FALSE, 
                       htmt2 = FALSE)
  htmt_cor <- calchtmt(data = data, model = model_est, latent1 = latent1, latent2 = latent2, scale = TRUE, 
                       htmt2 = FALSE)
  htmt_2_cor_value <- calchtmt(data = data, model = model_est, latent1 = latent1, latent2 = latent2, scale = TRUE, 
                               htmt2 = TRUE)
  
  return(list(mga = mga_value$out, 
              htmt_cov = htmt_cov$out, 
              htmt_cor = htmt_cor$out,
              htmt_2_cor = htmt_2_cor_value$out,
              mga_error = mga_value$error,
              mga_warning = mga_value$warning,
              ))
}


bs_p_value <- function(data, nboot, model, latent1, latent2, model_constrained, model_unconstrained)
{
  test_sta <- wrap(data = data, model = model, latent1 = latent1, latent2 = latent2)
  
  
  boot <- sapply(1:nboot, function(x) wrap(data = dplyr::sample_n(data, nrow(data), replace = TRUE), model = model, latent1 = latent1, latent2 = latent2))
  mga_boot <- unlist(boot["mga",])
  htmt_cov_boot <- unlist(boot["htmt_cov",])
  htmt_cor_boot <- unlist(boot["htmt_cor",])
  htmt_2_cor_boot <- unlist(boot["htmt_2_cor",])
  
  valid_mga_boot <- mga_boot[!is.na(mga_boot)]
  valid_htmt_cov_boot <- htmt_cov_boot[!is.na(htmt_cov_boot)]
  valid_htmt_cor_boot <- htmt_cor_boot[!is.na(htmt_cor_boot)]
  valid_htmt_2_cor_boot <- htmt_2_cor_boot[!is.na(htmt_2_cor_boot)]
  
  p_mga <- length(valid_mga_boot[valid_mga_boot > 1]) / length(valid_mga_boot)
  p_htmt_cov <- length(valid_htmt_cov_boot[valid_htmt_cov_boot > 1]) / length(valid_htmt_cov_boot)
  p_htmt_cor <- length(valid_htmt_cor_boot[valid_htmt_cor_boot > 1]) / length(valid_htmt_cor_boot)
  p_htmt_2 <- length(valid_htmt_2_cor_boot[valid_htmt_2_cor_boot > 1]) / length(valid_htmt_2_cor_boot)
  
  con_phi <- constrained_phi(model_constrained = model_constrained, model_unconstrained, data = data)
  
  FL <- fornell_larcker(data = data, model = model, latent1 = latent1, latent2 = latent2)

  
  return(list(test_stat = list(mga = test_sta$mga, 
                               htmt_cov = test_sta$htmt_cov,
                               htmt_cor = test_sta$htmt_cor,
                               htmt_2 = test_sta$htmt_2_cor,
                               con_phi = con_phi$test_stat), 
              p_values = list(p_mga = p_mga, 
                              p_htmt_cov = p_htmt_cov, 
                              p_htmt_cor = p_htmt_cor, 
                              p_htmt_2 = p_htmt_2, 
                              p_con_phi = con_phi$p), 
              boot = boot, 
              FL_dec =  FL$dec,
              FL_avesv = FL$ave_sv
  )) 
}


wrapper_function <- function(data, model, latent1, latent2, model_constrained, model_unconstrained){

  mga_value <- multigroup(data = data, model = model, 
                          latent1 = latent1, latent2 = latent2)
  htmt_cov <- calchtmt(data = data, model = model, latent1 = latent1, 
                       latent2 = latent2, scale = FALSE, htmt2 = FALSE)
  htmt_cor <- calchtmt(data = data, model = model, latent1 = latent1, 
                       latent2 = latent2, scale = TRUE,htmt2 = FALSE)
  htmt_2_cor_value <- calchtmt(data = data, model = model, latent1 = latent1, 
                               latent2 = latent2, scale = TRUE, htmt2 = TRUE)
 
  con_phi <- constrained_phi(data = data, 
                             model_constrained = model_constrained, 
                             model_unconstrained = model_unconstrained
                             )

  FL <- fornell_larcker(data, model, latent1, latent2)

    
  output <- data.frame(test = c("MGA", "HTMT_COV", "HTMT_COR", "HTMT_2", "CON_PHI", "FL"), 
                       stat = c(mga_value$out, htmt_cov$out, htmt_cor$out, htmt_2_cor_value$out, NA, NA),
                       constrained = c(NA, NA, NA, NA, con_phi$constrained, NA), 
                       unconstrained = c(NA, NA, NA, NA, con_phi$unconstrained, NA),
                       FL_1 = c(NA, NA, NA, NA, NA, FL$ave_sv[1]), 
                       FL_2 = c(NA, NA, NA, NA, NA, FL$ave_sv[2]),
                       warning = c( mga_value$warning, htmt_cov$warning, htmt_cor$warning, htmt_2_cor_value$warning, NA, FL$warning),
                       error =  c( mga_value$error, htmt_cov$error, htmt_cor$error, htmt_2_cor_value$error, NA, FL$error),
                       constrained_warning = c(NA, NA, NA, NA, con_phi$constrained_warning, NA),
                       unconstrained_warning = c(NA, NA, NA, NA, con_phi$unconstrained_warning, NA),
                       constrained_error = c(NA, NA, NA, NA, con_phi$constrained_error, NA),
                       unconstrained_error = c(NA, NA, NA, NA, con_phi$unconstrained_error, NA)
                       ) 
  return(output)
}


#bs_p_value(data, nboot = 500, model = model_est, latent1 = "xi_1", latent2 = "xi_2", model_unconstrained, model_constrained)





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
model_constrained <- '
              #  latent variables
                xi_1 =~ x11 + x12 + x13
                xi_2 =~ x21 + x22 + x23 
                
                xi_1 ~~ 1 * xi_2
                xi_1 ~~ 1 * xi_1
                xi_2 ~~ 1 * xi_2
              ' 

model_unconstrained <- '
              #  latent variables
                xi_1 =~ x11 + x12 + x13
                xi_2 =~ x21 + x22 + x23 
                
                xi_1 ~~ xi_2
              ' 
