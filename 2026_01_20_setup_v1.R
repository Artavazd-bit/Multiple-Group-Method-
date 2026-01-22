# MGA function
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
  
  Dinvsqrt <- diag(1/sqrt(diag(Lzero)))
  L <- Dinvsqrt %*% Lzero %*% Dinvsqrt
  Fmatrix <- t(Fzero) %*% Dinvsqrt
  return(L[2,1])
}

# nthroot  function
nthroot <- function(x, n){
  if(x <0 && n %% 2 == 1 ){
    sign(x) * abs(x)^(1/n)
  }else{
    x^(1/n)
  }
}

# HTMT function
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
  
    if(htmt2 == FALSE){
      A = 1/(K_i*K_j) * sum(cor_values$val[cor_values$type == "het"])
      B = 2/(K_i*(K_i-1)) *  sum(cor_values$val[cor_values$type == "mono1"]) 
      C = 2/(K_j*(K_j-1)) *  sum(cor_values$val[cor_values$type == "mono2"]) 
      HTMT <- A / ((B*C)^(1/2))
    }
    else if(htmt2 == TRUE){
      #A =  (prod(cor_values$val[cor_values$type == "het"]))^(1/(K_i*K_j))
      A = nthroot((prod(cor_values$val[cor_values$type == "het"])), (K_i*K_j))
      #B =  (prod(cor_values$val[cor_values$type == "mono1"]))^(2/(K_i*(K_i-1))) 
      B = nthroot((prod(cor_values$val[cor_values$type == "mono1"])), (K_i*(K_i-1))/2)
      #C =  (prod(cor_values$val[cor_values$type == "mono2"]))^(2/(K_j*(K_j-1))) 
      C = nthroot((prod(cor_values$val[cor_values$type == "mono2"])), (K_j*(K_j-1))/2 ) 
      HTMT <- A / ((B*C)^(1/2))
    }
    else{
      print("ERROR")
    }
  if(is.na(A)){
    warning("nominator is NaN")
  }
  if(B< 0){
    warning("monoblock1 is negative")
  }
  if(C< 0){
    warning("monoblock2 is negative")
  }
  return(HTMT)
}

#error handling function
jasonssuperdupererrorhandlingfunction <- function(fun, ...){
  warnings <- list()
  error <- NA
  
  result <- tryCatch(
    withCallingHandlers(
      fun(...),
      warning = function(cond){
        warnings <<- c(warnings, list(cond))
        invokeRestart("muffleWarning")
      }
    ), 
    error = function(cond){
     error <<- conditionMessage(cond)
     NA
    }
  )
  list(res = result, warnings = warnings, error = error)
}

# constrained phi function
c_phi_and_FL <- function(model_constrained, model_unconstrained, data, latent1, latent2){
  con_model <- jasonssuperdupererrorhandlingfunction(fun = sem, 
                                                     model = model_constrained, 
                                                     data = data
                                                     )
  uncon_model <- jasonssuperdupererrorhandlingfunction(fun = sem, 
                                                       model = model_unconstrained, 
                                                       data = data
                                                       )

  avesv <- AVE(uncon_model$res)[c(latent1, latent2)] / 
      lavInspect(uncon_model$res, "cor.lv")[latent1, latent2]
  
  
  return(list(
    con_ts = con_model$res@Fit@test$standard$stat,
    con_df = con_model$res@Fit@test$standard$df,
    
    uncon_ts = uncon_model$res@Fit@test$standard$stat,
    uncon_df = uncon_model$res@Fit@test$standard$df,
    
    con_warning = con_model$warnings,
    uncon_warning = uncon_model$warnings, 
    
    con_error = con_model$error,
    uncon_error = uncon_model$error,
    
    ave_cor_ratio_1 = avesv[1],
    ave_cor_ratio_2 = avesv[2] 
    #warning = out$warnings,
    #error = out$error
  ))
}

#Fornell Larcker function
FL <- function(data, model, latent1, latent2){
  out <- jasonssuperdupererrorhandlingfunction(fun = sem, model = model, data = data)
  ave <- AVE(out$res)
  ave_sel <- ave[c(latent1, latent2)]
  cormatrix <- lavInspect(out$res, "cor.lv")
  avesv <- ave_sel / cormatrix[latent1, latent2]
  return(list(
    ave1 = avesv[1],
    ave2 = avesv[2], 
    warning = out$warnings,
    error = out$error
  ))
}

run_methods <- function(data, model_unconstrained, latent1, latent2, model_constrained){
  mga <- jasonssuperdupererrorhandlingfunction(
    fun = multigroup, 
    data = data,
    model = model_unconstrained,
    latent1 = latent1,
    latent2 = latent2 
  )
  htmt_cov <- jasonssuperdupererrorhandlingfunction(
    fun = calchtmt, 
    data = data,
    model = model_unconstrained, 
    latent1 = latent1, 
    latent2 = latent2, 
    scale = FALSE, 
    htmt2 = FALSE
  )
  htmt_cor <- jasonssuperdupererrorhandlingfunction(
    fun = calchtmt, 
    data = data,
    model = model_unconstrained, 
    latent1 = latent1, 
    latent2 = latent2, 
    scale = TRUE, 
    htmt2 = FALSE
  )
  htmt_2_cor <- jasonssuperdupererrorhandlingfunction(
    fun = calchtmt, 
    data = data,
    model = model_unconstrained, 
    latent1 = latent1, 
    latent2 = latent2, 
    scale = TRUE, 
    htmt2 = TRUE
  )
  c_phi_and_FL_res <- c_phi_and_FL(model_constrained = model_constrained, 
                                   model_unconstrained = model_unconstrained, 
                                   data = data, 
                                   latent1 = latent1, 
                                   latent2 = latent2)
  
  output <- dplyr::bind_rows(
    data.frame(test = "mga", stat = mga$res, warning = I(list(mga$warnings)), error = I(list(mga$error)) ),
    data.frame(test = "htmt_cov", stat = htmt_cov$res, warning = I(list(htmt_cov$warnings)), error = I(list(htmt_cov$error))), 
    data.frame(test = "htmt_cor", stat = htmt_cor$res, warning = I(list(htmt_cor$warnings)), error = I(list(htmt_cor$error))),
    data.frame(test = "htmt_2", stat = htmt_2_cor$res, warning = I(list(htmt_2_cor$warnings)), error = I(list(htmt_2_cor$error))),
    
    data.frame(test = "conphi", chisq_constrained = c_phi_and_FL_res$con_ts, chisq_unconstrained = c_phi_and_FL_res$uncon_ts, 
               df_constrained = c_phi_and_FL_res$con_df, df_unconstrained = c_phi_and_FL_res$uncon_df, 
               warning_constrained = I(list(c_phi_and_FL_res$con_warning)), warning = I(list(c_phi_and_FL_res$uncon_warning)), 
               error_constrained = I(list(c_phi_and_FL_res$con_error)), error = I(list(c_phi_and_FL_res$uncon_error))),
    
    data.frame(test = "fl", ave_cor_ratio_1 = c_phi_and_FL_res$ave_cor_ratio_1, ave_cor_ratio_2 = c_phi_and_FL_res$ave_cor_ratio_1, 
               warning = I(list(c_phi_and_FL_res$uncon_warning)), error = I(list(c_phi_and_FL_res$uncon_error)))
  )
  return(output)
}

#sim models
{
  
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
}

model_constrained <- '
              #  latent variables
                xi_1 =~ NA*x11 + x12 + x13
                xi_2 =~ NA*x21 + x22 + x23 
                
                xi_1 ~~ 1 * xi_2
                xi_1 ~~ 1 * xi_1
                xi_2 ~~ 1 * xi_2
              ' 

model_unconstrained <- '
              #  latent variables
                xi_1 =~ NA*x11 + x12 + x13
                xi_2 =~ NA*x21 + x22 + x23 
                
                xi_1 ~~ xi_2
                xi_1 ~~ 1 * xi_1
                xi_2 ~~ 1 * xi_2
              ' 



