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
  error <- NULL
  
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


